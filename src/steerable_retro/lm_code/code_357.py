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
    Detects a synthesis strategy involving a benzothiophene-containing compound
    that remains intact throughout the synthesis.
    """
    # Initialize flags
    benzothiophene_count = 0
    has_consistent_benzothiophene = False

    def dfs_traverse(node):
        nonlocal benzothiophene_count

        if node["type"] == "mol":
            # Check for benzothiophene in molecules
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                benzothiophene_pattern = Chem.MolFromSmarts("c1ccc2c(c1)ccs2")
                if mol.HasSubstructMatch(benzothiophene_pattern):
                    benzothiophene_count += 1
                    print(f"Found benzothiophene in molecule: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # Strategy is present if benzothiophene is found in multiple molecules
    has_consistent_benzothiophene = benzothiophene_count >= 3

    if has_consistent_benzothiophene:
        print(
            f"Detected benzothiophene-containing synthesis strategy with {benzothiophene_count} occurrences"
        )
    else:
        print(
            f"Benzothiophene-containing synthesis strategy not detected, found only {benzothiophene_count} occurrences"
        )

    return has_consistent_benzothiophene
