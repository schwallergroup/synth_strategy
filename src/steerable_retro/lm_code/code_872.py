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
    This function detects a strategy where a quinoline core is progressively
    elaborated through sequential functionalization steps.
    """
    # Track if we found a quinoline core
    found_quinoline = False
    # Track functionalization steps
    functionalization_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal found_quinoline, functionalization_steps

        if node["type"] == "mol":
            # Check if molecule contains quinoline core
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                quinoline_pattern = Chem.MolFromSmarts(
                    "[#6]1[#6][#6][#6][#6]2[#7][#6][#6][#6][#6]21"
                )
                if mol.HasSubstructMatch(quinoline_pattern):
                    found_quinoline = True
                    print(f"Found quinoline core at depth {depth}")

        elif node["type"] == "reaction" and found_quinoline:
            # Count functionalization steps on the quinoline core
            functionalization_steps += 1
            print(f"Found functionalization step at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have a quinoline core with multiple functionalization steps
    if found_quinoline and functionalization_steps >= 4:
        print("Detected quinoline elaboration strategy with multiple functionalization steps")
        return True
    return False
