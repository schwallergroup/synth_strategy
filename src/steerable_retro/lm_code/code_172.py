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
    Detects if the synthetic route creates a biaryl system containing methoxy groups.
    """
    methoxy_biaryl_detected = False

    def dfs_traverse(node):
        nonlocal methoxy_biaryl_detected

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for methoxy group
                methoxy_pattern = Chem.MolFromSmarts("c[O][CH3]")
                has_methoxy = mol.HasSubstructMatch(methoxy_pattern)

                # Check for biaryl system
                biaryl_pattern = Chem.MolFromSmarts("c:c-c:c")
                has_biaryl = mol.HasSubstructMatch(biaryl_pattern)

                if has_methoxy and has_biaryl:
                    print("Detected methoxy-containing biaryl system")
                    methoxy_biaryl_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return methoxy_biaryl_detected
