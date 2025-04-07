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

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    Detects if the final step (depth 0) is a Boc deprotection
    """
    boc_deprotection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_deprotection_found

        print(f"Examining node at depth {depth}, type: {node['type']}")

        # Check if this is a reaction node
        if node["type"] == "reaction":
            print(f"Found reaction node at depth {depth}")

            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Reaction SMILES: {rsmi}")

                # Check if this is a Boc deprotection reaction using the checker
                is_boc_deprotection = checker.check_reaction("Boc amine deprotection", rsmi)

                if not is_boc_deprotection:
                    # Try other Boc deprotection reaction types
                    is_boc_deprotection = checker.check_reaction(
                        "Boc amine deprotection of guanidine", rsmi
                    )

                if not is_boc_deprotection:
                    is_boc_deprotection = checker.check_reaction(
                        "Boc amine deprotection to NH-NH2", rsmi
                    )

                if not is_boc_deprotection:
                    is_boc_deprotection = checker.check_reaction(
                        "Tert-butyl deprotection of amine", rsmi
                    )

                if is_boc_deprotection:
                    print(f"Found Boc deprotection reaction at depth {depth}")
                    # If this is the final reaction (depth 0 or 1)
                    if depth <= 1:
                        boc_deprotection_found = True
                        print("This is a late-stage Boc deprotection")
                else:
                    print("Not a standard Boc deprotection reaction, checking manually")

                    # Additional verification by checking for Boc group disappearance
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        print(f"Reactants: {reactants}")
                        print(f"Product: {product}")

                        # Check if any reactant has a Boc group that's not in the product
                        boc_in_reactants = False
                        for reactant in reactants:
                            has_boc = checker.check_fg("Boc", reactant)
                            print(f"Reactant {reactant} has Boc: {has_boc}")
                            if has_boc:
                                boc_in_reactants = True
                                break

                        boc_in_product = checker.check_fg("Boc", product)
                        print(f"Product has Boc: {boc_in_product}")

                        if boc_in_reactants and not boc_in_product:
                            print("Found Boc group removal (manual check)")
                            # If this is the final reaction (depth 0 or 1)
                            if depth <= 1:
                                boc_deprotection_found = True
                                print("This is a late-stage Boc deprotection")
                    except Exception as e:
                        print(f"Error during manual Boc check: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting traversal of synthesis route")
    dfs_traverse(route)
    print(f"Final result: {boc_deprotection_found}")

    return boc_deprotection_found
