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
    This function detects a strategy that is primarily linear but ends with a convergent step
    (typically a cross-coupling reaction).
    """
    # Track reaction depths and types
    reaction_types = {}  # depth -> reaction type
    max_depth = [0]  # Use a list to allow modification in nested function

    def calculate_depth(node, parent_depth=0):
        """Calculate the depth of a node in the synthesis tree"""
        if node["type"] == "reaction":
            depth = parent_depth
            max_depth[0] = max(max_depth[0], depth)
            return depth
        else:  # It's a molecule node
            if not node.get("children", []):  # Leaf node (starting material)
                return parent_depth
            max_child_depth = parent_depth
            for child in node.get("children", []):
                child_depth = calculate_depth(child, parent_depth + 1)
                max_child_depth = max(max_child_depth, child_depth)
            return max_child_depth

    # First pass: calculate depths
    calculate_depth(route)

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Count number of significant reactants (ignoring small molecules)
            significant_reactants = 0
            for r in reactants_smiles:
                mol = Chem.MolFromSmiles(r)
                if mol and mol.GetNumHeavyAtoms() > 5:  # Consider only substantial fragments
                    significant_reactants += 1

            # Check if this is a cross-coupling reaction
            is_cross_coupling = (
                checker.check_reaction("Suzuki", rsmi)
                or checker.check_reaction("Negishi", rsmi)
                or checker.check_reaction("Stille", rsmi)
                or checker.check_reaction("Heck", rsmi)
                or checker.check_reaction("Sonogashira", rsmi)
                or checker.check_reaction("Buchwald-Hartwig", rsmi)
                or checker.check_reaction("Ullmann", rsmi)
            )

            # Classify the reaction
            if is_cross_coupling and significant_reactants >= 2:
                reaction_type = "convergent"
            elif significant_reactants >= 2:
                reaction_type = "convergent"
            else:
                reaction_type = "linear"

            reaction_types[depth] = reaction_type
            print(
                f"Reaction at depth {depth} classified as {reaction_type} with {significant_reactants} significant reactants"
            )

            # Continue traversal with incremented depth
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)
        else:  # It's a molecule node
            for child in node.get("children", []):
                dfs_traverse(child, depth)

    # Second pass: classify reactions
    dfs_traverse(route)

    # Check if the pattern matches linear-to-convergent
    if not reaction_types:
        print("No reactions found in the route")
        return False

    # The final step (depth 0) should be convergent
    final_step_convergent = reaction_types.get(0) == "convergent"

    # Most earlier steps should be linear
    linear_steps = sum(
        1 for depth, type in reaction_types.items() if depth > 0 and type == "linear"
    )
    total_earlier_steps = sum(1 for depth in reaction_types if depth > 0)

    # Strategy is true if final step is convergent and at least 70% of earlier steps are linear
    linear_ratio = linear_steps / max(1, total_earlier_steps) if total_earlier_steps > 0 else 0
    is_linear_to_convergent = final_step_convergent and (linear_ratio >= 0.7)

    print(f"Final step convergent: {final_step_convergent}")
    print(f"Linear steps: {linear_steps}/{total_earlier_steps} ({linear_ratio:.2f})")
    print(f"Linear-to-convergent strategy: {is_linear_to_convergent}")

    return is_linear_to_convergent
