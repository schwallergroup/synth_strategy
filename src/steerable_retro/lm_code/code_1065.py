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
    This function detects if the synthetic route follows a convergent pattern
    where multiple fragments are prepared separately and then joined.
    """
    # Track branches in the synthesis
    branches = []

    def count_branches(node, branch_id=0):
        """Count the number of separate branches in the synthesis"""
        if node["type"] == "reaction":
            reaction_smiles = node["metadata"]["rsmi"]
            reactants = reaction_smiles.split(">")[0].split(".")

            # If there are multiple reactants, this might be a convergent step
            if len(reactants) > 1:
                # For each child, assign a new branch ID
                new_branches = []
                for i, child in enumerate(node.get("children", [])):
                    new_branch_id = branch_id * 10 + (i + 1)  # Create unique branch IDs
                    child_branches = count_branches(child, new_branch_id)
                    new_branches.extend(child_branches)
                return new_branches
            else:
                # Continue with same branch ID
                for child in node.get("children", []):
                    return count_branches(child, branch_id)

        # If this is a molecule with no children, it's a leaf node (starting material)
        if node["type"] == "mol" and not node.get("children"):
            return [branch_id]

        # Process children
        all_branches = []
        for child in node.get("children", []):
            child_branches = count_branches(child, branch_id)
            all_branches.extend(child_branches)

        return all_branches

    # Start branch counting
    branches = count_branches(route)

    # If we have multiple unique branches, it's a convergent synthesis
    unique_branches = set(branches)
    is_convergent = len(unique_branches) > 1

    if is_convergent:
        print(f"Found convergent synthesis with {len(unique_branches)} branches")

    return is_convergent
