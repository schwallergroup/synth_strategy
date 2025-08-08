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
    This function detects if the synthesis preserves a benzothiazole core throughout.
    """
    benzothiazole_count = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal benzothiazole_count, total_reactions

        if node["type"] == "reaction":
            total_reactions += 1
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check if product contains benzothiazole
                benzothiazole_pattern = Chem.MolFromSmarts(
                    "[#6]1[#6][#6][#6][#6]2[#6]1[#7][#6][#16]2"
                )

                product_mol = Chem.MolFromSmiles(product)
                if product_mol is not None and product_mol.HasSubstructMatch(benzothiazole_pattern):
                    benzothiazole_count += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Check if benzothiazole is present in all or most reactions
    preserved = (benzothiazole_count > 0) and (benzothiazole_count >= total_reactions * 0.8)
    if preserved:
        print(
            f"Detected benzothiazole preservation throughout synthesis ({benzothiazole_count}/{total_reactions} reactions)"
        )

    return preserved
