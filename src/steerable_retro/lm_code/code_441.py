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
    This function detects if a tertiary alcohol group is preserved throughout the synthesis.
    """
    has_tertiary_alcohol_in_final = False
    has_tertiary_alcohol_in_intermediate = False

    def dfs_traverse(node, depth=0):
        nonlocal has_tertiary_alcohol_in_final, has_tertiary_alcohol_in_intermediate

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # SMARTS pattern for tertiary alcohol
                tert_alcohol_pattern = Chem.MolFromSmarts("[C]([C])([C])[OH]")
                if mol.HasSubstructMatch(tert_alcohol_pattern):
                    if depth == 0:
                        has_tertiary_alcohol_in_final = True
                        print(f"Found tertiary alcohol in final product: {node['smiles']}")
                    else:
                        has_tertiary_alcohol_in_intermediate = True
                        print(
                            f"Found tertiary alcohol in intermediate at depth {depth}: {node['smiles']}"
                        )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_tertiary_alcohol_in_final and has_tertiary_alcohol_in_intermediate
