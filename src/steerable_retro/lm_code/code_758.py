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
    This function detects if key structural motifs (difluorophenyl and methylpyridine-cyano)
    are maintained throughout the synthesis.
    """
    # Track presence of motifs at different depths
    difluorophenyl_depths = set()
    methylpyridine_cyano_depths = set()

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for difluorophenyl motif
                difluorophenyl_pattern = Chem.MolFromSmarts("c1c(F)c(F)ccc1")
                if mol.HasSubstructMatch(difluorophenyl_pattern):
                    difluorophenyl_depths.add(depth)
                    print(f"Difluorophenyl motif found at depth {depth}")

                # Check for methylpyridine-cyano motif
                methylpyridine_cyano_pattern = Chem.MolFromSmarts("[CH3]c1cc(C#N)cnc1")
                if mol.HasSubstructMatch(methylpyridine_cyano_pattern):
                    methylpyridine_cyano_depths.add(depth)
                    print(f"Methylpyridine-cyano motif found at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if both motifs are present at multiple depths
    difluorophenyl_maintained = len(difluorophenyl_depths) > 1
    methylpyridine_cyano_maintained = len(methylpyridine_cyano_depths) > 1

    result = difluorophenyl_maintained and methylpyridine_cyano_maintained
    print(f"Maintained structural motifs strategy: {result}")
    print(f"Difluorophenyl found at depths: {difluorophenyl_depths}")
    print(f"Methylpyridine-cyano found at depths: {methylpyridine_cyano_depths}")

    return result
