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
    This function detects late-stage isoxazole formation strategy.
    It looks for the creation of an isoxazole ring in the final or penultimate synthetic step.
    """
    isoxazole_formed = False

    def dfs_traverse(node, depth=0):
        nonlocal isoxazole_formed

        # Check final step (depth 0) and penultimate step (depth 1)
        if node["type"] == "reaction" and depth <= 1:
            # Extract reactants and product
            try:
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if not rsmi:
                    print(f"No reaction SMILES found in metadata at depth {depth}")
                    return

                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")
                print(f"Product: {product_smiles}")
                print(f"Reactants: {reactants_smiles}")

                # Check if isoxazole is in product
                if checker.check_ring("isoxazole", product_smiles):
                    print(f"Found isoxazole ring in product at depth {depth}")

                    # Check if isoxazole is not in any reactant
                    isoxazole_in_reactants = False
                    for reactant in reactants_smiles:
                        if checker.check_ring("isoxazole", reactant):
                            print(f"Isoxazole found in reactant: {reactant}")
                            isoxazole_in_reactants = True
                            break

                    if not isoxazole_in_reactants:
                        print(
                            f"Isoxazole not found in reactants at depth {depth} - confirming late-stage formation"
                        )

                        # Check if this is a known isoxazole formation reaction
                        isoxazole_reaction = False

                        # List of reactions that can form isoxazoles
                        isoxazole_forming_reactions = [
                            "Huisgen 1,3 dipolar cycloaddition",
                            "[3+2]-cycloaddition of hydrazone and alkyne",
                            "[3+2]-cycloaddition of hydrazone and alkene",
                            "[3+2]-cycloaddition of diazoalkane and alkyne",
                            "[3+2]-cycloaddition of diazoalkane and alkene",
                            "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                            "Huisgen alkene-azide 1,3 dipolar cycloaddition",
                            "oxadiazole",
                            "pyrazole formation",
                        ]

                        for reaction_type in isoxazole_forming_reactions:
                            if checker.check_reaction(reaction_type, rsmi):
                                print(
                                    f"Confirmed isoxazole formation reaction: {reaction_type}"
                                )
                                isoxazole_reaction = True
                                break

                        # Even if not a known reaction type, if isoxazole appears in product but not reactants,
                        # we'll consider it an isoxazole formation
                        if isoxazole_reaction or depth <= 1:
                            print(
                                f"Isoxazole formed at depth {depth} - this is late-stage formation"
                            )
                            isoxazole_formed = True
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    print(f"Final result: isoxazole_formed = {isoxazole_formed}")

    return isoxazole_formed
