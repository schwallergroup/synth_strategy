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
    This function detects a linear fragment assembly strategy where three or more
    distinct fragments are combined sequentially.
    """
    # Track the number of fragment combinations
    fragment_combinations = 0

    def dfs_traverse(node):
        nonlocal fragment_combinations

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # If there are multiple reactants, this is a fragment combination
            if len(reactants_smiles) >= 2:
                fragment_combinations += 1
                print(f"Found fragment combination with {len(reactants_smiles)} reactants")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return true if we found at least 3 fragment combinations
    if fragment_combinations >= 2:  # At least 3 fragments total (2 combinations)
        print(
            f"Linear fragment assembly strategy detected with {fragment_combinations} combinations"
        )
        return True

    return False
