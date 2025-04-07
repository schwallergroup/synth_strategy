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


def main(route):
    """
    Detects if the synthesis follows a primarily linear rather than convergent approach.
    This is determined by analyzing the structure of the synthetic tree.
    """
    # Track branching in the tree
    max_branches = 0
    convergent_pattern_found = False

    def is_significant_reactant(smiles):
        """Determine if a reactant is significant (not just a reagent)"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False

            # Check molecular complexity
            num_atoms = mol.GetNumAtoms()
            num_rings = len(Chem.GetSSSR(mol))

            # Simple reagents typically have few atoms and no rings
            if num_atoms <= 5 and num_rings == 0:
                return False

            # Check for common reagents
            common_reagents = ["O", "C(=O)O", "CC(=O)O", "Cl", "Br", "I", "F", "[OH-]", "[H+]"]
            if smiles in common_reagents:
                return False

            return True
        except:
            # If there's an error, assume it's significant
            return True

    def dfs_traverse(node, depth=0, branch_path=None):
        nonlocal max_branches, convergent_pattern_found

        if branch_path is None:
            branch_path = []

        if node["type"] == "reaction":
            # Count number of significant reactants in this reaction
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                significant_reactants = [r for r in reactants if is_significant_reactant(r)]
                num_significant_reactants = len(significant_reactants)

                print(
                    f"Depth {depth}: Found {num_significant_reactants} significant reactants out of {len(reactants)} total"
                )

                # Update max branches if needed
                if num_significant_reactants > max_branches:
                    max_branches = num_significant_reactants

                # Check for convergent pattern - multiple complex reactants coming together
                if num_significant_reactants > 2:
                    convergent_pattern_found = True
                    print(
                        f"Convergent pattern found at depth {depth} with {num_significant_reactants} significant reactants"
                    )

        # Process children
        child_count = len(node.get("children", []))

        # If this node has multiple children, it might indicate branching
        if child_count > 1 and node["type"] == "mol":
            print(f"Branching detected at depth {depth} with {child_count} paths")

        for i, child in enumerate(node.get("children", [])):
            new_branch_path = branch_path + [i]
            dfs_traverse(child, depth + 1, new_branch_path)

    # Start traversal
    dfs_traverse(route)

    # Linear synthesis has max 2 significant reactants per step and no convergent patterns
    is_linear = max_branches <= 2 and not convergent_pattern_found
    print(
        f"Max significant branches: {max_branches}, Convergent pattern: {convergent_pattern_found}"
    )
    return is_linear
