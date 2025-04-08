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
    This function detects if the synthesis involves formation of sulfonamide
    functional groups from amines.
    """
    sulfonamide_formation_found = False

    def dfs_traverse(node):
        nonlocal sulfonamide_formation_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if the reaction is a known sulfonamide formation reaction
            if checker.check_reaction(
                "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
            ) or checker.check_reaction(
                "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
            ):
                print(f"Detected sulfonamide formation via Schotten-Baumann reaction: {rsmi}")
                sulfonamide_formation_found = True
                return

            # If not a known reaction type, check for the transformation manually
            product_mol = Chem.MolFromSmiles(product)

            # Check if product contains sulfonamide
            if product_mol and checker.check_fg("Sulfonamide", product):
                print(f"Product contains sulfonamide: {product}")

                # Check if any reactant contains sulfonyl chloride
                sulfonyl_chloride_present = False
                for reactant in reactants:
                    if checker.check_fg("Sulfonyl halide", reactant):
                        print(f"Reactant contains sulfonyl halide: {reactant}")
                        sulfonyl_chloride_present = True
                        break

                # Check if any reactant contains amine (primary or secondary)
                amine_present = False
                for reactant in reactants:
                    if (
                        checker.check_fg("Primary amine", reactant)
                        or checker.check_fg("Secondary amine", reactant)
                        or checker.check_fg("Aniline", reactant)
                    ):
                        print(f"Reactant contains amine: {reactant}")
                        amine_present = True
                        break

                # Verify that the sulfonamide is formed in this reaction
                if sulfonyl_chloride_present and amine_present:
                    # Check if sulfonamide is not present in any reactant
                    sulfonamide_in_reactants = False
                    for reactant in reactants:
                        if checker.check_fg("Sulfonamide", reactant):
                            sulfonamide_in_reactants = True
                            break

                    if not sulfonamide_in_reactants:
                        print(f"Sulfonamide formation detected in reaction: {rsmi}")
                        sulfonamide_formation_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Sulfonamide formation strategy detected: {sulfonamide_formation_found}")
    return sulfonamide_formation_found
