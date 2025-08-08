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
    Detects a strategy where a specific functional group (hydroxymethyl)
    is preserved throughout the synthesis while the molecule is being built.
    """
    hydroxymethyl_pattern = Chem.MolFromSmarts("[CH2][OH]")
    depths_with_hydroxymethyl = set()
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(hydroxymethyl_pattern):
                depths_with_hydroxymethyl.add(depth)
                print(f"Hydroxymethyl group found at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if hydroxymethyl group is present at multiple depths
    if len(depths_with_hydroxymethyl) <= 1:
        return False

    # Check if hydroxymethyl is preserved from early to late stages
    has_early = any(d > max_depth / 2 for d in depths_with_hydroxymethyl)
    has_late = any(d <= max_depth / 2 for d in depths_with_hydroxymethyl)

    preserved_throughout = has_early and has_late

    print(f"Preserved hydroxymethyl group strategy detected: {preserved_throughout}")
    print(f"Depths with hydroxymethyl: {depths_with_hydroxymethyl}")

    return preserved_throughout
