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
    This function detects if a cyano group is preserved throughout the synthesis.
    """
    cyano_present_all_steps = True

    def dfs_traverse(node):
        nonlocal cyano_present_all_steps

        if node["type"] == "mol" and node["smiles"]:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                cyano_pattern = Chem.MolFromSmarts("[C]#[N]")

                if mol and not mol.HasSubstructMatch(cyano_pattern):
                    # If any molecule doesn't have cyano group, set flag to False
                    if not node.get("in_stock", False):  # Ignore starting materials
                        cyano_present_all_steps = False
                        print(f"Molecule without cyano group found: {node['smiles']}")
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return cyano_present_all_steps
