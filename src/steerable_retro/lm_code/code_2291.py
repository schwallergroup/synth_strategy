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
    This function detects thiazole ring formation from alpha-bromo ketone and thioacetamide.
    """
    thiazole_formation_detected = False

    def dfs_traverse(node):
        nonlocal thiazole_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")

                # Check if the reaction is a thiazole formation reaction
                is_thiazole_reaction = checker.check_reaction("thiazole", rsmi)

                # Check if product contains thiazole ring
                has_thiazole_in_product = checker.check_ring("thiazole", product)

                if has_thiazole_in_product:
                    print(f"Product contains thiazole ring: {product}")

                # Check if reactants contain alpha-bromo ketone and thioamide
                has_alpha_bromo_ketone = False
                has_thioamide = False

                for reactant in reactants:
                    # Check for alpha-bromo ketone (ketone + primary or secondary halide)
                    if checker.check_fg("Ketone", reactant) and (
                        checker.check_fg("Primary halide", reactant)
                        or checker.check_fg("Secondary halide", reactant)
                    ):
                        has_alpha_bromo_ketone = True
                        print(f"Found alpha-bromo ketone: {reactant}")

                    # Check for alpha-bromo ketone pattern in SMILES
                    if "C(=O)CH2Br" in reactant or "[CH2:21]Br" in reactant:
                        has_alpha_bromo_ketone = True
                        print(f"Found alpha-bromo ketone pattern: {reactant}")

                    # Check for thioamide
                    if checker.check_fg("Thioamide", reactant):
                        has_thioamide = True
                        print(f"Found thioamide: {reactant}")

                # If we can't directly check for thioamide, look for thiourea or similar structures
                if not has_thioamide:
                    for reactant in reactants:
                        if (
                            checker.check_fg("Thiourea", reactant)
                            or "C(=S)N" in reactant
                            or "[C:2]([NH2:3])=[S:22]" in reactant
                        ):
                            has_thioamide = True
                            print(f"Found thioamide-like structure: {reactant}")

                # Determine if this is a thiazole formation reaction
                if (
                    is_thiazole_reaction
                    or (has_thiazole_in_product and has_alpha_bromo_ketone and has_thioamide)
                    or (
                        has_thiazole_in_product
                        and ("CH2Br" in rsmi or "[CH2:21]Br" in rsmi)
                        and ("C(=S)N" in rsmi or "[C:2]([NH2:3])=[S:22]" in rsmi)
                    )
                ):
                    print(
                        f"Detected thiazole formation from alpha-bromo ketone and thioacetamide\nReactants: {reactants}\nProduct: {product}"
                    )
                    thiazole_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return thiazole_formation_detected
