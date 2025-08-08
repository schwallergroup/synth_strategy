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
    Detects if a methoxy group on an aromatic ring is preserved throughout the synthesis.
    """
    # Track methoxy aromatic systems through the synthesis
    methoxy_systems = []
    preservation_confirmed = False

    def dfs_traverse(node, depth=0):
        nonlocal methoxy_systems, preservation_confirmed

        if node["type"] == "mol":
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                # Check for methoxy on aromatic pattern
                methoxy_pattern = Chem.MolFromSmarts("[#6]:[#6][#8][#6]")
                if mol.HasSubstructMatch(methoxy_pattern):
                    methoxy_systems.append((depth, smiles))

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if methoxy system appears at multiple depths
    depth_counts = len(set([d for d, _ in methoxy_systems]))

    if depth_counts >= 3:  # Present in at least 3 different depths
        preservation_confirmed = True
        print(f"Methoxy group preserved across {depth_counts} synthesis steps")
        for depth, smiles in methoxy_systems:
            print(f"  Depth {depth}: {smiles}")

    return preservation_confirmed
