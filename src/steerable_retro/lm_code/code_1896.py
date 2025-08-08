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
    This function detects a linear synthesis strategy where a biaryl system
    (difluorophenyl-isoxazole) is preserved throughout the synthesis.
    """
    # Track if we consistently see the biaryl system
    biaryl_present_count = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal biaryl_present_count, total_reactions

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]
            total_reactions += 1

            # Check for difluorophenyl-isoxazole biaryl system
            difluorophenyl_pattern = Chem.MolFromSmarts("[F]c1[cH][c]([F])[cH][cH]c1")
            isoxazole_pattern = Chem.MolFromSmarts("c1conc1")

            product_mol = Chem.MolFromSmiles(product) if product else None

            if (
                product_mol
                and product_mol.HasSubstructMatch(difluorophenyl_pattern)
                and product_mol.HasSubstructMatch(isoxazole_pattern)
            ):
                biaryl_present_count += 1
                print(f"Detected difluorophenyl-isoxazole biaryl system in: {product}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if the biaryl system is preserved throughout (in at least 80% of reactions)
    result = total_reactions > 0 and (biaryl_present_count / total_reactions) >= 0.8
    print(f"Linear synthesis with preserved biaryl strategy detected: {result}")
    print(f"Biaryl present in {biaryl_present_count} out of {total_reactions} reactions")
    return result
