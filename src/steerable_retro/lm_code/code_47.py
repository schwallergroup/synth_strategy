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
    Detects if the synthetic route follows a linear synthesis strategy
    rather than a convergent approach.
    """
    # Track the maximum branching factor in the synthesis tree
    max_branching = 0

    def dfs_traverse(node):
        nonlocal max_branching

        if node["type"] == "reaction":
            # Count the number of molecular children (reactants)
            mol_children = [child for child in node.get("children", []) if child["type"] == "mol"]
            branching = len(mol_children)
            max_branching = max(max_branching, branching)

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If max branching is <= 2, it's likely a linear synthesis
    # (one main reactant + one reagent/catalyst)
    is_linear = max_branching <= 2
    if is_linear:
        print("Detected linear synthesis strategy")
    else:
        print(f"Detected convergent synthesis strategy with branching factor {max_branching}")

    return is_linear
