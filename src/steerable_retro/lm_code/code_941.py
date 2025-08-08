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
    This function detects late-stage formation of morpholine amides in the synthetic route.
    Late stage is defined as occurring in the first half of the synthesis (lower depth values).
    """
    morpholine_amide_pattern = Chem.MolFromSmarts("[#6]-C(=O)-N1CCOCC1")
    max_depth = 0
    morpholine_amide_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal max_depth, morpholine_amide_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(morpholine_amide_pattern):
                    if morpholine_amide_depth is None or depth < morpholine_amide_depth:
                        morpholine_amide_depth = depth
                        print(f"Morpholine amide formation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if morpholine amide formation occurs in the first half of synthesis
    if morpholine_amide_depth is not None and morpholine_amide_depth <= max_depth / 2:
        return True
    return False
