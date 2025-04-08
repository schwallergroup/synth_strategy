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
    Detects if the synthesis involves formation of a triazolone ring from an acyl chloride.
    """
    triazolone_formed = False
    acyl_chloride_present = False

    def dfs_traverse(node):
        nonlocal triazolone_formed, acyl_chloride_present

        if node["type"] == "mol":
            # Check if molecule contains triazolone
            if node["smiles"]:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    triazolone_pattern = Chem.MolFromSmarts("[#7]1[#7][#7][#7][#6]1=[#8]")
                    if mol.HasSubstructMatch(triazolone_pattern):
                        triazolone_formed = True

                    # Check for acyl chloride
                    acyl_chloride_pattern = Chem.MolFromSmarts("[#6](=[#8])[Cl]")
                    if mol.HasSubstructMatch(acyl_chloride_pattern):
                        acyl_chloride_present = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Triazolone formation from acyl chloride: {triazolone_formed and acyl_chloride_present}")
    return triazolone_formed and acyl_chloride_present
