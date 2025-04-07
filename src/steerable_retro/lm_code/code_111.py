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
    This function detects if the synthesis follows a linear strategy with sequential transformations.
    It checks if most reactions have only one main reactant (1:1 transformations).
    """
    reaction_count = 0
    linear_reaction_count = 0

    # Track the depth of each reaction to analyze tree structure
    max_depth = 0
    reaction_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, linear_reaction_count, max_depth

        # Update max depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Count significant reactants (ignoring very small molecules)
            significant_reactants = 0
            for reactant_smiles in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant_smiles)
                if (
                    mol and mol.GetNumHeavyAtoms() > 3
                ):  # Consider molecules with >3 heavy atoms significant
                    significant_reactants += 1

            reaction_count += 1
            reaction_depths.append(depth)

            # Check if this is a linear transformation (exactly one significant reactant)
            if significant_reactants == 1:
                linear_reaction_count += 1

            # Check if this is a known linear reaction type
            elif (
                checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rsmi)
                or checker.check_reaction("Reduction of aldehydes and ketones to alcohols", rsmi)
                or checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                or checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                )
            ):
                linear_reaction_count += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Calculate linearity metrics
    if reaction_count > 0:
        # Reaction linearity: percentage of reactions that are linear transformations
        reaction_linearity = linear_reaction_count / reaction_count

        # Path linearity: how close the synthesis is to a single path
        # In a perfectly linear synthesis, max_depth + 1 = reaction_count
        path_linearity = (max_depth + 1) / reaction_count if reaction_count > 0 else 0

        # Combined linearity score (weighted average)
        linearity_score = 0.7 * reaction_linearity + 0.3 * path_linearity

        print(
            f"Linear reactions: {linear_reaction_count}/{reaction_count} ({reaction_linearity:.2f})"
        )
        print(
            f"Path linearity: {path_linearity:.2f} (max depth: {max_depth}, reactions: {reaction_count})"
        )
        print(f"Overall linearity score: {linearity_score:.2f}")

        # Consider it linear if the combined score is at least 0.65
        return linearity_score >= 0.65

    return False
