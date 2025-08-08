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
    This function detects a linear synthesis strategy with sequential transformations
    (as opposed to convergent synthesis with multiple fragments).
    """
    reaction_count = 0
    max_children_per_reaction = 0

    def dfs_traverse(node):
        nonlocal reaction_count, max_children_per_reaction

        if node["type"] == "reaction":
            reaction_count += 1
            # Count number of reactant children
            reactant_count = sum(1 for child in node.get("children", []) if child["type"] == "mol")
            max_children_per_reaction = max(max_children_per_reaction, reactant_count)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # A linear synthesis typically has at most 2 reactants per reaction
    is_linear = reaction_count > 0 and max_children_per_reaction <= 2

    print(
        f"Linear synthesis strategy detected: {is_linear} (max reactants per step: {max_children_per_reaction})"
    )
    return is_linear
