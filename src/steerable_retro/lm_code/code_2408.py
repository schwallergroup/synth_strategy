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
    Detects if the synthesis route involves a late-stage aromatic halogenation
    (addition of halogen to aromatic ring in the final or penultimate step).
    """
    halogenation_detected = False

    def get_depth(node, current_depth=0):
        """Calculate the depth of a node in the synthesis route"""
        if node == route:  # Root node
            return 0

        # Traverse the tree to find the node and its depth
        def find_node_depth(current_node, target, depth=0):
            if current_node == target:
                return depth

            for child in current_node.get("children", []):
                result = find_node_depth(child, target, depth + 1)
                if result is not None:
                    return result
            return None

        return find_node_depth(route, node, 0)

    def dfs_traverse(node, depth=0):
        nonlocal halogenation_detected

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            # Check if this is a late-stage reaction (depth 0 or 1)
            if depth <= 1:
                rsmi = node["metadata"]["rsmi"]

                # Check if this is an aromatic halogenation reaction
                is_halogenation = (
                    checker.check_reaction("Aromatic fluorination", rsmi)
                    or checker.check_reaction("Aromatic chlorination", rsmi)
                    or checker.check_reaction("Aromatic bromination", rsmi)
                    or checker.check_reaction("Aromatic iodination", rsmi)
                )

                if is_halogenation:
                    print(f"Late-stage aromatic halogenation detected at depth {depth}")
                    halogenation_detected = True
                else:
                    # Alternative check: verify if halogen is added to aromatic ring
                    try:
                        reactants_part = rsmi.split(">")[0]
                        product_part = rsmi.split(">")[-1]

                        # Check if product has aromatic halide but reactants don't
                        reactants_have_aromatic_halide = any(
                            checker.check_fg("Aromatic halide", r)
                            for r in reactants_part.split(".")
                        )

                        product_has_aromatic_halide = checker.check_fg(
                            "Aromatic halide", product_part
                        )

                        if product_has_aromatic_halide and not reactants_have_aromatic_halide:
                            print(
                                f"Late-stage aromatic halogenation detected at depth {depth} (FG check)"
                            )
                            halogenation_detected = True
                    except Exception as e:
                        print(f"Error analyzing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return halogenation_detected
