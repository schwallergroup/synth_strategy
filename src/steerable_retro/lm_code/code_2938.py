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
    This function detects if the synthesis uses a convergent approach in the final step,
    combining 3 or more fragments.
    """
    final_step_reactant_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal final_step_reactant_count

        # Check if this is a reaction node
        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            # Check if this is the final step (depth 0)
            if depth == 0:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                # Count non-empty reactants
                valid_reactants = [r for r in reactants if r.strip()]
                final_step_reactant_count = len(valid_reactants)
                print(f"Final step has {final_step_reactant_count} reactants: {valid_reactants}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root (target molecule)
    if route["type"] == "mol":
        # If the root is a molecule, check its children for the final reaction step
        for child in route.get("children", []):
            dfs_traverse(child)
    else:
        # If the root is already a reaction, start from there
        dfs_traverse(route)

    print(f"Final reactant count: {final_step_reactant_count}")
    return final_step_reactant_count >= 3
