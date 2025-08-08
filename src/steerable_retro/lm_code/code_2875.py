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
    This function detects a strategy involving epoxide formation followed by epoxide opening.
    """
    epoxide_formation = False
    epoxide_opening = False

    def dfs_traverse(node):
        nonlocal epoxide_formation, epoxide_opening

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for epoxide formation
            epoxide_pattern = Chem.MolFromSmarts("[C]1[O][C]1")

            # Check if epoxide is in product but not in reactants
            if product and product.HasSubstructMatch(epoxide_pattern):
                reactants_have_epoxide = any(
                    r and r.HasSubstructMatch(epoxide_pattern) for r in reactants if r
                )
                if not reactants_have_epoxide:
                    epoxide_formation = True
                    print("Detected epoxide formation")

            # Check for epoxide opening
            if any(r and r.HasSubstructMatch(epoxide_pattern) for r in reactants if r):
                if not (product and product.HasSubstructMatch(epoxide_pattern)):
                    epoxide_opening = True
                    print("Detected epoxide opening")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Return True if both epoxide formation and opening are detected
    return epoxide_formation and epoxide_opening
