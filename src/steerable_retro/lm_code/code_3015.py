#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    This function detects a linear synthesis pattern (as opposed to convergent).

    A linear synthesis has a single main transformation path where each step
    produces one significant intermediate that becomes the reactant for the next step.
    """
    # Track the structure of the synthesis tree
    branch_counts = []

    def dfs_traverse(node, depth=0):
        # Extend branch_counts list if needed
        while len(branch_counts) <= depth:
            branch_counts.append(0)

        if node["type"] == "mol":
            # Skip leaf nodes (starting materials)
            if node.get("in_stock", False):
                return

            # Count significant branches at this depth
            significant_children = 0
            for child in node.get("children", []):
                if child["type"] == "reaction":
                    # For each reaction, count significant reactants
                    if "metadata" in child and "rsmi" in child["metadata"]:
                        rsmi = child["metadata"]["rsmi"]
                        reactants_smiles = rsmi.split(">")[0].split(".")

                        significant_reactants = 0
                        for r in reactants_smiles:
                            mol = Chem.MolFromSmiles(r)
                            if mol and mol.GetNumHeavyAtoms() > 3:  # Significant reactant
                                significant_reactants += 1

                        if significant_reactants > 0:
                            significant_children += 1

            # Update branch count at this depth
            branch_counts[depth] = max(branch_counts[depth], significant_children)

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Filter out zeros (depths with no branching)
    non_zero_branches = [count for count in branch_counts if count > 0]

    # A linear synthesis should have no more than one significant branch at any depth
    is_linear = all(count <= 1 for count in non_zero_branches)

    # Check if we have enough data to make a determination
    if not non_zero_branches:
        print("Warning: Could not determine synthesis pattern - insufficient data")
        return False

    if is_linear:
        print(f"Detected linear synthesis pattern with branch counts: {non_zero_branches}")
    else:
        print(f"Detected convergent synthesis pattern with branch counts: {non_zero_branches}")

    return is_linear
