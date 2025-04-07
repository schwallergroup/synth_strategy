#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    Detects if the synthesis involves formation of a triazole ring through
    a hydrazide cyclization sequence.
    """
    has_triazole_formation = False

    def dfs_traverse(node):
        nonlocal has_triazole_formation

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            # Check for hydrazide pattern in reactants
            hydrazide_pattern = Chem.MolFromSmarts("[NH][NH][C]=[O]")

            # Check for triazole pattern in products
            triazole_pattern = Chem.MolFromSmarts("[n]1[c][nH][n][c]1")

            try:
                reactants_mol = Chem.MolFromSmiles(reactants_part)
                has_hydrazide = reactants_mol and reactants_mol.HasSubstructMatch(
                    hydrazide_pattern
                )
            except:
                has_hydrazide = False

            try:
                product_mol = Chem.MolFromSmiles(products_part)
                has_triazole = product_mol and product_mol.HasSubstructMatch(
                    triazole_pattern
                )
            except:
                has_triazole = False

            if has_hydrazide and has_triazole:
                print(
                    f"Detected triazole formation via hydrazide cyclization at depth {node['metadata'].get('depth', 'unknown')}"
                )
                has_triazole_formation = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_triazole_formation
