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
    This function detects a functional group interconversion sequence:
    ester → carboxylic acid → acyl chloride → ketone
    """
    # SMARTS patterns for functional groups
    ester_pattern = Chem.MolFromSmarts("[#6][CX3](=[OX1])[OX2][#6]")
    acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H]")
    acyl_chloride_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[Cl]")
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=[OX1])[#6]")

    # Track if we found each transformation
    found_ester = False
    found_acid = False
    found_acyl_chloride = False
    found_ketone = False

    def dfs_traverse(node):
        nonlocal found_ester, found_acid, found_acyl_chloride, found_ketone

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                if mol.HasSubstructMatch(ester_pattern):
                    found_ester = True
                    print("Found ester functional group")
                if mol.HasSubstructMatch(acid_pattern):
                    found_acid = True
                    print("Found carboxylic acid functional group")
                if mol.HasSubstructMatch(acyl_chloride_pattern):
                    found_acyl_chloride = True
                    print("Found acyl chloride functional group")
                if mol.HasSubstructMatch(ketone_pattern):
                    found_ketone = True
                    print("Found ketone functional group")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we found the complete sequence
    result = found_ester and found_acid and found_acyl_chloride and found_ketone
    print(f"Functional group interconversion sequence: {result}")
    return result
