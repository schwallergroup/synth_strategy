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
    This function detects a synthetic strategy involving multiple S-N and C-N bond
    manipulations in a linear sequence.
    """
    # Initialize counters
    sn_bond_changes = 0
    cn_bond_changes = 0

    def dfs_traverse(node):
        nonlocal sn_bond_changes, cn_bond_changes

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check for S-N bond changes
            if re.search(r"\[S.*\].*\[N", rsmi) or re.search(r"\[N.*\].*\[S", rsmi):
                sn_bond_changes += 1
                print(f"Found S-N bond manipulation in reaction: {rsmi}")

            # Check for C-N bond changes
            if re.search(r"\[C.*\].*\[N", rsmi) or re.search(r"\[N.*\].*\[C", rsmi):
                cn_bond_changes += 1
                print(f"Found C-N bond manipulation in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we have at least one S-N and one C-N bond manipulation
    strategy_present = sn_bond_changes >= 1 and cn_bond_changes >= 1

    print(f"S-N and C-N bond manipulation strategy detected: {strategy_present}")
    print(f"S-N bond manipulations: {sn_bond_changes}")
    print(f"C-N bond manipulations: {cn_bond_changes}")

    return strategy_present
