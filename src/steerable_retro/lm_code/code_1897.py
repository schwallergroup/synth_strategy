#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

root_data = "/home/andres/Documents/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    This function detects thiazole ring formation from acyclic precursors.
    """
    thiazole_formation_detected = False

    def dfs_traverse(node):
        nonlocal thiazole_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if thiazole is in product
                product_has_thiazole = checker.check_ring("thiazole", product_smiles)

                if product_has_thiazole:
                    print(f"Product has thiazole: {product_smiles}")

                    # Check if thiazole is in any reactant
                    thiazole_in_reactants = False
                    for reactant in reactants_smiles:
                        if checker.check_ring("thiazole", reactant):
                            thiazole_in_reactants = True
                            print(f"Thiazole found in reactant: {reactant}")
                            break

                    # If thiazole is in product but not in reactants, it's a thiazole formation
                    if not thiazole_in_reactants:
                        print(f"Potential thiazole formation detected in reaction: {rsmi}")

                        # Check if reactants are acyclic (don't contain thiazole or related rings)
                        acyclic_precursors = True
                        for reactant in reactants_smiles:
                            if (
                                checker.check_ring("thiazole", reactant)
                                or checker.check_ring("benzothiazole", reactant)
                                or checker.check_ring("isothiazole", reactant)
                            ):
                                acyclic_precursors = False
                                print(f"Reactant contains thiazole-related ring: {reactant}")
                                break

                        if acyclic_precursors:
                            # Check for all possible thiazole formation reactions
                            if (
                                checker.check_reaction("{thiazole}", rsmi)
                                or checker.check_reaction("thiazole", rsmi)
                                or checker.check_reaction("benzothiazole", rsmi)
                                or checker.check_reaction("{benzothiazole}", rsmi)
                                or checker.check_reaction(
                                    "benzothiazole_derivatives_aldehyde", rsmi
                                )
                                or checker.check_reaction(
                                    "benzothiazole_derivatives_carboxylic-acid/ester", rsmi
                                )
                                or checker.check_reaction("Formation of NOS Heterocycles", rsmi)
                            ):
                                print(f"Thiazole formation reaction confirmed: {rsmi}")
                                thiazole_formation_detected = True
                            else:
                                # Check if it's a thiazole formation even if not explicitly categorized
                                # Look for common functional groups in reactants that could form thiazole
                                has_thiol_or_thioamide = False
                                has_amino_or_amide = False
                                has_carbonyl = False

                                for reactant in reactants_smiles:
                                    if (
                                        checker.check_fg("Aliphatic thiol", reactant)
                                        or checker.check_fg("Aromatic thiol", reactant)
                                        or checker.check_fg("Thioamide", reactant)
                                        or checker.check_fg("Thiourea", reactant)
                                    ):
                                        has_thiol_or_thioamide = True

                                    if (
                                        checker.check_fg("Primary amine", reactant)
                                        or checker.check_fg("Secondary amine", reactant)
                                        or checker.check_fg("Primary amide", reactant)
                                        or checker.check_fg("Secondary amide", reactant)
                                    ):
                                        has_amino_or_amide = True

                                    if (
                                        checker.check_fg("Aldehyde", reactant)
                                        or checker.check_fg("Ketone", reactant)
                                        or checker.check_fg("Carboxylic acid", reactant)
                                        or checker.check_fg("Ester", reactant)
                                        or checker.check_fg("Acyl halide", reactant)
                                    ):
                                        has_carbonyl = True

                                if has_thiol_or_thioamide and (has_amino_or_amide or has_carbonyl):
                                    print(
                                        f"Thiazole formation detected based on functional groups: {rsmi}"
                                    )
                                    thiazole_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return thiazole_formation_detected
