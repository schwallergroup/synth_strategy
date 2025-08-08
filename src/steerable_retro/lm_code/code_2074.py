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
    Detects a synthetic route that includes ketone reduction to an alcohol.
    """
    has_ketone_reduction = False

    def dfs_traverse(node):
        nonlocal has_ketone_reduction

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for ketone reduction
            ketone_pattern = Chem.MolFromSmarts("C(=O)")
            alcohol_pattern = Chem.MolFromSmarts("C([OH])")

            if any(
                mol.HasSubstructMatch(ketone_pattern) for mol in reactants
            ) and product.HasSubstructMatch(alcohol_pattern):
                print("Found ketone reduction")
                has_ketone_reduction = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Ketone reduction strategy detected: {has_ketone_reduction}")
    return has_ketone_reduction
