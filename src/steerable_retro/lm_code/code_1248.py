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
    Detects if the synthesis uses a late-stage amide reduction strategy.
    Specifically looks for an amide reduction in the final step (depth 0).
    """
    print("Starting late_stage_amide_reduction_strategy analysis")

    # Validate route structure
    if not isinstance(route, dict) or "type" not in route:
        print("Invalid route structure")
        return False

    # Ensure the root node is a molecule (the target compound)
    if route["type"] != "mol":
        print(f"Root node is not a molecule, but {route['type']}")
        return False

    amide_reduction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_reduction_found

        print(f"Traversing node of type {node['type']} at depth {depth}")

        if node["type"] == "mol":
            print(f"Molecule node: {node.get('smiles', 'No SMILES')}")

        # Check if this is a reaction node
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Reaction node at depth {depth}: {rsmi}")

            # Check if this is the final step (depth 0)
            if depth == 1:  # The reaction leading to the target molecule
                print(f"Checking final step reaction: {rsmi}")

                # Check if this is an amide reduction reaction using the checker functions
                if (
                    checker.check_reaction("Reduction of primary amides to amines", rsmi)
                    or checker.check_reaction("Reduction of secondary amides to amines", rsmi)
                    or checker.check_reaction("Reduction of tertiary amides to amines", rsmi)
                ):
                    print(f"Found amide reduction reaction in final step: {rsmi}")
                    amide_reduction_found = True
                else:
                    # If specific reaction check fails, try a more general approach
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        print(f"Reactants: {reactants}")
                        print(f"Product: {product}")

                        # Check for amide in reactants and amine in product
                        amide_in_reactants = False
                        for reactant in reactants:
                            if (
                                checker.check_fg("Primary amide", reactant)
                                or checker.check_fg("Secondary amide", reactant)
                                or checker.check_fg("Tertiary amide", reactant)
                            ):
                                amide_in_reactants = True
                                print(f"Found amide in reactant: {reactant}")
                                break

                        amine_in_product = (
                            checker.check_fg("Primary amine", product)
                            or checker.check_fg("Secondary amine", product)
                            or checker.check_fg("Tertiary amine", product)
                        )

                        if amine_in_product:
                            print(f"Found amine in product: {product}")

                        # If we have both amide in reactants and amine in product, it's likely an amide reduction
                        if amide_in_reactants and amine_in_product:
                            print("Detected amide reduction pattern in final step")
                            amide_reduction_found = True
                    except Exception as e:
                        print(f"Error analyzing reaction: {e}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root node (target molecule)
    dfs_traverse(route)
    print(f"Amide reduction found: {amide_reduction_found}")
    return amide_reduction_found
