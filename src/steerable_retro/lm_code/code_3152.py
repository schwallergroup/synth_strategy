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
    Detects if brominated aromatic core is preserved throughout the synthesis.
    """
    brominated_aromatic_count = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal brominated_aromatic_count, total_reactions

        if node["type"] == "reaction":
            total_reactions += 1

            # Extract product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            product_smiles = parts[2]

            # Check if product contains brominated aromatic
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None
            if product:
                brominated_aromatic_pattern = Chem.MolFromSmarts("c[Br]")
                if product.HasSubstructMatch(brominated_aromatic_pattern):
                    brominated_aromatic_count += 1
                    print(f"Found brominated aromatic in product")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if brominated aromatic is preserved in all reactions
    is_preserved = brominated_aromatic_count == total_reactions and total_reactions > 0
    print(f"Brominated aromatic preservation detected: {is_preserved}")
    return is_preserved
