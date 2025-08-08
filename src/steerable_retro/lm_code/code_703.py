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
    This function detects if the synthesis follows a linear strategy
    where each reaction builds upon a single product from the previous step.
    """
    is_linear = True
    reaction_count = 0
    max_reactants_per_reaction = 0

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, reaction_count, max_reactants_per_reaction

        if node["type"] == "reaction":
            reaction_count += 1
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                num_reactants = len(reactants)
                max_reactants_per_reaction = max(max_reactants_per_reaction, num_reactants)

                # If any reaction has more than 2 reactants, it's likely not a linear synthesis
                if num_reactants > 2:
                    is_linear = False
                    print(
                        f"Found non-linear reaction with {num_reactants} reactants at depth {depth}"
                    )

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # A linear synthesis should have multiple reactions and mostly 1-2 reactants per reaction
    if reaction_count < 2:
        is_linear = False

    print(
        f"Synthesis has {reaction_count} reactions with max {max_reactants_per_reaction} reactants per reaction"
    )
    return is_linear
