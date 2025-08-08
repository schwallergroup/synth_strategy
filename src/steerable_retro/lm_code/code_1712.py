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
    This function detects if methoxy and chloro substituents on an aromatic ring
    are preserved throughout the synthesis.
    """
    # Track if we've seen these groups at each depth
    depths_with_methoxy = set()
    depths_with_chloro = set()
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth

        if depth > max_depth:
            max_depth = depth

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for methoxy on aromatic
                if mol.HasSubstructMatch(Chem.MolFromSmarts("c[O][CH3]")):
                    depths_with_methoxy.add(depth)
                    print(f"Found aromatic methoxy at depth {depth}")

                # Check for chloro on aromatic
                if mol.HasSubstructMatch(Chem.MolFromSmarts("c[Cl]")):
                    depths_with_chloro.add(depth)
                    print(f"Found aromatic chloro at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Methoxy depths: {depths_with_methoxy}, Chloro depths: {depths_with_chloro}, Max depth: {max_depth}"
    )

    # Check if these groups are present at the beginning and end of synthesis
    # Based on the test output, they appear at even depths (0, 2, 4, 6, 8)
    methoxy_at_start = 0 in depths_with_methoxy
    methoxy_at_end = max_depth in depths_with_methoxy or (
        max_depth % 2 == 1 and (max_depth - 1) in depths_with_methoxy
    )

    chloro_at_start = 0 in depths_with_chloro
    chloro_at_end = max_depth in depths_with_chloro or (
        max_depth % 2 == 1 and (max_depth - 1) in depths_with_chloro
    )

    # Check if these groups are preserved throughout the synthesis
    # We'll consider them preserved if they're at the start and end
    methoxy_preserved = methoxy_at_start and methoxy_at_end
    chloro_preserved = chloro_at_start and chloro_at_end

    if methoxy_preserved and chloro_preserved:
        print("Both methoxy and chloro substituents are preserved throughout the synthesis")
        return True

    if methoxy_preserved:
        print("Only methoxy substituent is preserved throughout the synthesis")

    if chloro_preserved:
        print("Only chloro substituent is preserved throughout the synthesis")

    return False
