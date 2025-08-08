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
    This function detects a linear synthesis strategy with sequential addition of fragments
    rather than convergent synthesis.
    """
    reaction_count = 0
    max_reactants_per_step = 0

    def dfs_traverse(node):
        nonlocal reaction_count, max_reactants_per_step

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if rsmi:
                reactants_smiles = rsmi.split(">")[0].split(".")
                num_reactants = len(reactants_smiles)
                max_reactants_per_step = max(max_reactants_per_step, num_reactants)
                reaction_count += 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Linear synthesis typically has 2 or fewer reactants per step
    # and at least 2 reaction steps
    is_linear = max_reactants_per_step <= 2 and reaction_count >= 2

    if is_linear:
        print(
            f"Detected linear synthesis with {reaction_count} steps and max {max_reactants_per_step} reactants per step"
        )

    return is_linear
