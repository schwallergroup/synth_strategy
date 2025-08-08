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
    This function detects a linear synthesis strategy with a late-stage fragment coupling.
    """
    # Track the number of fragments combined at each depth
    fragment_couplings = {}
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")

            # Count number of fragments being combined
            num_fragments = len(reactants_smiles)
            if num_fragments > 1:
                fragment_couplings[depth] = num_fragments
                print(f"Found fragment coupling at depth {depth} with {num_fragments} fragments")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the strategy is present:
    # 1. Linear synthesis: most steps have only one fragment
    # 2. Late-stage coupling: fragment coupling at depth 0 or 1
    has_late_coupling = 0 in fragment_couplings or 1 in fragment_couplings
    is_mostly_linear = (
        len(fragment_couplings) <= max_depth / 2
    )  # Less than half of steps are couplings

    strategy_present = has_late_coupling and is_mostly_linear

    if strategy_present:
        print("Detected linear synthesis with late-stage fragment coupling")

    return strategy_present
