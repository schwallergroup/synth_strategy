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
    Detects if the synthesis follows a linear pattern without convergent steps.
    """
    max_reactants_per_step = 0

    def dfs_traverse(node):
        nonlocal max_reactants_per_step

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Count number of reactants
            num_reactants = len(reactants_smiles)
            max_reactants_per_step = max(max_reactants_per_step, num_reactants)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Consider it linear if no step has more than 2 reactants
    is_linear = max_reactants_per_step <= 2

    if is_linear:
        print(
            f"Confirmed linear synthesis pattern (max reactants per step: {max_reactants_per_step})"
        )
    else:
        print(f"Detected convergent synthesis (max reactants per step: {max_reactants_per_step})")

    return is_linear
