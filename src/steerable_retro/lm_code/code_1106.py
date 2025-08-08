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
    Detects if the synthesis maintains the core scaffold while modifying functional groups.
    This is characterized by preserving aromatic systems throughout the synthesis.
    """
    aromatic_systems = []
    scaffold_preserved = True

    def dfs_traverse(node):
        nonlocal aromatic_systems, scaffold_preserved

        if node["type"] == "mol":
            smiles = node["smiles"]
            if smiles:
                # Count aromatic rings
                aromatic_count = smiles.count("c")
                if aromatic_count > 0:
                    if not aromatic_systems:
                        aromatic_systems.append(aromatic_count)
                    else:
                        # If the number of aromatic carbons changes significantly,
                        # the scaffold is not preserved
                        if abs(aromatic_count - aromatic_systems[0]) > 2:
                            print(
                                f"Scaffold changed: aromatic count from {aromatic_systems[0]} to {aromatic_count}"
                            )
                            scaffold_preserved = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return scaffold_preserved and len(aromatic_systems) > 0
