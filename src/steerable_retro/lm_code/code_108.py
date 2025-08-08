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
    This function detects if the synthesis follows a linear strategy (each step builds on previous).

    A linear synthesis is defined as a route where:
    1. Each reaction has exactly one non-stock molecule that continues the synthesis path
    2. There are multiple reaction steps
    """
    # Track the number of reaction steps and if we have a non-linear path
    reaction_steps = 0
    is_non_linear = False

    def dfs_traverse(node, depth=0):
        nonlocal reaction_steps, is_non_linear

        if node["type"] == "reaction":
            reaction_steps += 1
            print(f"Found reaction at depth {depth}")

            # Get all molecule children
            mol_children = [child for child in node.get("children", []) if child["type"] == "mol"]

            # Count non-stock molecules that continue the synthesis
            continuing_mols = []
            for mol_child in mol_children:
                # A molecule continues synthesis if it has reaction children
                has_reaction_children = any(
                    grandchild["type"] == "reaction" for grandchild in mol_child.get("children", [])
                )
                if has_reaction_children and not mol_child.get("in_stock", False):
                    continuing_mols.append(mol_child)

            # In a linear synthesis, at most one non-stock molecule should continue
            if len(continuing_mols) > 1:
                is_non_linear = True
                print(
                    f"Non-linear path detected at depth {depth} with {len(continuing_mols)} continuing molecules"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root (target molecule)
    dfs_traverse(route)

    print(f"Total reaction steps: {reaction_steps}")

    # A linear synthesis has multiple steps and is not non-linear
    is_linear = reaction_steps > 1 and not is_non_linear

    if is_linear:
        print("Detected linear synthesis strategy")
    else:
        print("Not a linear synthesis strategy")

    return is_linear
