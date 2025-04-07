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
    This function detects if the synthesis involves manipulation of a
    thiophene-containing heterocyclic system.
    """
    thiophene_present = False

    def dfs_traverse(node):
        nonlocal thiophene_present

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    thiophene_pattern = Chem.MolFromSmarts("c1ccsc1")
                    if mol.HasSubstructMatch(thiophene_pattern):
                        thiophene_present = True
                        print(f"Thiophene detected in molecule: {node['smiles']}")
            except:
                pass
        elif node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]
            try:
                mol = Chem.MolFromSmiles(product_smiles)
                if mol:
                    thiophene_pattern = Chem.MolFromSmarts("c1ccsc1")
                    if mol.HasSubstructMatch(thiophene_pattern):
                        thiophene_present = True
                        print(
                            f"Thiophene detected in reaction product: {product_smiles}"
                        )
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return thiophene_present
