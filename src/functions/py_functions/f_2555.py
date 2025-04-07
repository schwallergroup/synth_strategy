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
    This function detects if the synthetic route involves compounds with multiple methoxy groups.
    """
    methoxy_rich = False

    def dfs_traverse(node):
        nonlocal methoxy_rich

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                methoxy_pattern = Chem.MolFromSmarts("[c]-[O][CH3]")
                if mol:
                    matches = mol.GetSubstructMatches(methoxy_pattern)
                    if len(matches) >= 3:  # At least 3 methoxy groups
                        print(
                            f"Methoxy-rich structure detected with {len(matches)} methoxy groups"
                        )
                        methoxy_rich = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return methoxy_rich
