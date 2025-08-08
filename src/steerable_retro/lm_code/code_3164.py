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
    Detects a synthetic strategy involving pyridine to piperidine transformation.
    """
    pyridine_found = False
    piperidine_found = False

    def dfs_traverse(node):
        nonlocal pyridine_found, piperidine_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for pyridine in reactants
            reactant_mol = Chem.MolFromSmiles(reactants)
            if reactant_mol and reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[n;r6]1ccccc1")):
                pyridine_found = True
                print("Found pyridine pattern in reactants")

            # Check for piperidine in products
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[N;r6]1CCCCC1")):
                piperidine_found = True
                print("Found piperidine pattern in products")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if both patterns were found in the route
    result = pyridine_found and piperidine_found
    print(f"Pyridine to piperidine strategy detected: {result}")
    return result
