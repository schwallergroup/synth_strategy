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
    This function detects if Boc protection is maintained throughout the synthesis.
    """
    # Track Boc presence in each reaction
    reactions_with_boc = []
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal reactions_with_boc, total_reactions

        if node["type"] == "reaction":
            total_reactions += 1

            # Extract product
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for Boc group
            boc_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]([C])([C])[C]")
            if product is not None and product.HasSubstructMatch(boc_pattern):
                reactions_with_boc.append(True)
                print(f"Detected Boc group in reaction product: {product_smiles}")
            else:
                reactions_with_boc.append(False)
                print(f"No Boc group in reaction product: {product_smiles}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if Boc is present in all reactions
    all_have_boc = all(reactions_with_boc) if reactions_with_boc else False

    print(f"Reactions with Boc: {sum(reactions_with_boc)}/{total_reactions}")
    print(f"All reactions have Boc: {all_have_boc}")

    return all_have_boc and total_reactions > 0
