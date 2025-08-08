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
    Detects if the synthesis uses a convergent approach with multiple fragments
    being combined rather than a purely linear synthesis.
    """
    fragment_combinations = 0

    def dfs_traverse(node):
        nonlocal fragment_combinations

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            reactants = rsmi.split(">")[0].split(".")

            # If a reaction has multiple reactants, it's combining fragments
            if len(reactants) > 1:
                print(f"Found fragment combination with {len(reactants)} reactants")
                fragment_combinations += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return fragment_combinations >= 2  # At least 2 fragment combinations for convergent synthesis
