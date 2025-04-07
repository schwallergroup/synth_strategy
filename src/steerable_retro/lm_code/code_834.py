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
    This function detects if the synthetic route involves a sequence of
    functional group transformations: aldehyde → alcohol → ester → nitrile
    """
    # SMARTS patterns for functional groups
    aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
    alcohol_pattern = Chem.MolFromSmarts("[CH2][OH]")
    ester_pattern = Chem.MolFromSmarts("[#6][C](=O)[O][#6]")
    nitrile_pattern = Chem.MolFromSmarts("[C]#N")

    # Track the sequence of functional groups observed
    sequence = []

    def dfs_traverse(node):
        nonlocal sequence

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for functional groups
                if mol.HasSubstructMatch(aldehyde_pattern) and "aldehyde" not in sequence:
                    sequence.append("aldehyde")
                    print(f"Aldehyde detected: {node['smiles']}")
                elif mol.HasSubstructMatch(alcohol_pattern) and "alcohol" not in sequence:
                    sequence.append("alcohol")
                    print(f"Alcohol detected: {node['smiles']}")
                elif mol.HasSubstructMatch(ester_pattern) and "ester" not in sequence:
                    sequence.append("ester")
                    print(f"Ester detected: {node['smiles']}")
                elif mol.HasSubstructMatch(nitrile_pattern) and "nitrile" not in sequence:
                    sequence.append("nitrile")
                    print(f"Nitrile detected: {node['smiles']}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if the sequence contains all four functional groups in the correct order
    correct_sequence = ["aldehyde", "alcohol", "ester", "nitrile"]

    # Check if all elements are present
    all_present = all(fg in sequence for fg in correct_sequence)

    # Check if the order is preserved (allowing for other elements in between)
    order_preserved = True
    last_idx = -1
    for fg in correct_sequence:
        if fg in sequence:
            idx = sequence.index(fg)
            if idx <= last_idx:
                order_preserved = False
                break
            last_idx = idx

    result = all_present and order_preserved
    print(f"Functional group sequence detected: {sequence}")
    print(f"Correct sequence found: {result}")

    return result
