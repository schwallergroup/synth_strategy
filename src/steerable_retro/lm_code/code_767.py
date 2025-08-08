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
    Detects if the route uses orthogonal amine protection strategies
    (Boc and phthalimide).
    """
    has_boc = False
    has_phthalimide = False

    def dfs_traverse(node, depth=0):
        nonlocal has_boc, has_phthalimide

        if node["type"] == "reaction":
            # Extract product
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]
            product = Chem.MolFromSmiles(product_smiles)

            # Define patterns
            boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)[N]")
            phthalimide_pattern = Chem.MolFromSmarts("O=C1c2ccccc2C(=O)N1")

            # Check for Boc and phthalimide
            if product and product.HasSubstructMatch(boc_pattern):
                print(f"Boc protection detected at depth {depth}")
                has_boc = True

            if product and product.HasSubstructMatch(phthalimide_pattern):
                print(f"Phthalimide protection detected at depth {depth}")
                has_phthalimide = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if both protection strategies are used
    return has_boc and has_phthalimide
