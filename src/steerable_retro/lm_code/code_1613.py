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
    This function detects if the synthesis follows a linear strategy where each intermediate
    is used exactly once in the next reaction.
    """
    # Track all molecules and their usage count
    molecule_usage = {}

    def dfs_traverse(node):
        if node["type"] == "mol":
            smiles = node["smiles"]
            if smiles in molecule_usage:
                molecule_usage[smiles] += 1
            else:
                molecule_usage[smiles] = 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if all intermediates are used exactly once (except starting materials and final product)
    is_linear = all(
        count == 1
        for smiles, count in molecule_usage.items()
        if not (
            smiles == route["smiles"]
            or any(
                node.get("in_stock", False)
                for node in route.get("children", [])
                if node["type"] == "mol" and node["smiles"] == smiles
            )
        )
    )

    if is_linear:
        print("Detected linear synthesis strategy: each intermediate used exactly once")

    return is_linear
