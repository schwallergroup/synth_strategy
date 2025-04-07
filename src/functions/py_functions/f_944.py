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
    Detects the presence of fluorinated aromatic rings in the synthetic route.
    """
    fluorinated_aromatics_count = 0

    def dfs_traverse(node):
        nonlocal fluorinated_aromatics_count

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for fluorinated aromatic pattern
                pattern = Chem.MolFromSmarts("c-[F]")
                matches = mol.GetSubstructMatches(pattern)
                if matches:
                    fluorinated_aromatics_count += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    if fluorinated_aromatics_count > 0:
        print(
            f"Fluorinated aromatics detected ({fluorinated_aromatics_count} instances)"
        )
        return True
    return False
