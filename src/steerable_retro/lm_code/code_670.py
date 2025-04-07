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
    Detects if the synthetic route involves a Boc deprotection as the final step.
    """
    print("Starting Boc deprotection strategy analysis...")

    if not route or "type" not in route:
        print("Invalid route structure")
        return False

    final_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal final_deprotection

        print(f"Analyzing node at depth {depth}, type: {node['type']}")

        # Check if this is a reaction node at the first step (depth 1)
        if node["type"] == "reaction" and depth == 1 and len(node.get("children", [])) > 0:
            print(f"Examining final reaction node at depth {depth}")

            try:
                if "metadata" in node and "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]
                    print(f"Reaction SMILES: {rsmi}")

                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Reactants: {reactants}")
                    print(f"Product: {product}")

                    # Check if this is a Boc deprotection reaction using various reaction types
                    boc_reaction_types = [
                        "Boc amine deprotection",
                        "Boc amine deprotection of guanidine",
                        "Boc amine deprotection to NH-NH2",
                        "Tert-butyl deprotection of amine",
                    ]

                    for rxn_type in boc_reaction_types:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Detected {rxn_type} in final step: {rsmi}")
                            final_deprotection = True
                            return

                    # If reaction check failed, verify by checking functional groups
                    if not final_deprotection:
                        print("Checking for Boc group in reactants and product...")
                        # Check if any reactant has Boc group but product doesn't
                        has_boc_in_reactants = any(checker.check_fg("Boc", r) for r in reactants)
                        has_boc_in_product = checker.check_fg("Boc", product)

                        # Check if product has amine group (expected after Boc removal)
                        has_amine_in_product = checker.check_fg(
                            "Primary amine", product
                        ) or checker.check_fg("Secondary amine", product)

                        print(
                            f"Boc in reactants: {has_boc_in_reactants}, Boc in product: {has_boc_in_product}"
                        )
                        print(f"Amine in product: {has_amine_in_product}")

                        if has_boc_in_reactants and not has_boc_in_product and has_amine_in_product:
                            print(f"Detected Boc group removal in final step: {rsmi}")
                            final_deprotection = True
                else:
                    print("No reaction SMILES found in metadata")
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {final_deprotection}")
    return final_deprotection
