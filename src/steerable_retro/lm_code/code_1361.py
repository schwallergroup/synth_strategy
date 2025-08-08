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
    This function detects if the synthetic route involves heterocyclic structures
    like benzimidazole and pyridine.
    """
    heterocycles_found = False

    def dfs_traverse(node):
        nonlocal heterocycles_found

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for benzimidazole
                benzimidazole_pattern = Chem.MolFromSmarts("c1nc2ccccc2[nH]1")
                # Check for pyridine
                pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")

                if mol.HasSubstructMatch(benzimidazole_pattern) or mol.HasSubstructMatch(
                    pyridine_pattern
                ):
                    print("Heterocyclic structure detected")
                    heterocycles_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return heterocycles_found
