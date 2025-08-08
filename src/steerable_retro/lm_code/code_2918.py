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
    Detects a synthetic strategy involving isocyanate as a key building block
    in the early stages of the synthesis.
    """
    has_isocyanate = False

    def dfs_traverse(node):
        nonlocal has_isocyanate

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for isocyanate in reactants or products
            isocyanate_pattern = Chem.MolFromSmarts("[#6][N]=[C]=[O]")

            for smiles in reactants_smiles + [product_smiles]:
                mol = Chem.MolFromSmiles(smiles)
                if mol and mol.HasSubstructMatch(isocyanate_pattern):
                    print(f"Detected isocyanate in molecule: {smiles}")
                    has_isocyanate = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_isocyanate
