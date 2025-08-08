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
    This function detects a synthetic strategy that builds complexity around a pyridine scaffold.
    """
    pyridine_core = False
    elaboration_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal pyridine_core, elaboration_steps

        if node["type"] == "mol":
            # Check if molecule contains pyridine core
            if "smiles" in node:
                try:
                    mol = Chem.MolFromSmiles(node["smiles"])
                    if mol:
                        pyridine_pattern = Chem.MolFromSmarts(
                            "[n;r6]1[c;r6][c;r6][c;r6][c;r6][c;r6]1"
                        )
                        if mol.HasSubstructMatch(pyridine_pattern):
                            pyridine_core = True
                except:
                    print("Error processing SMILES in pyridine detection")

        elif node["type"] == "reaction":
            if pyridine_core:
                elaboration_steps += 1

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Return True if we have a pyridine core and at least 2 elaboration steps
    result = pyridine_core and elaboration_steps >= 2
    if result:
        print(f"Pyridine elaboration strategy detected with {elaboration_steps} steps")

    return result
