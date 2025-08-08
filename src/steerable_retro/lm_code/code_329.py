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
    This function detects a linear fragment assembly strategy where multiple fragments
    are sequentially coupled without convergent steps.
    """
    # Track the number of fragment couplings
    fragment_couplings = 0
    # Track if the synthesis is linear (no convergent steps)
    is_linear = True

    def dfs_traverse(node):
        nonlocal fragment_couplings, is_linear

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # If there are multiple reactants, it's a fragment coupling
            if len(reactants_smiles) > 1:
                fragment_couplings += 1

                # Check if any of the reactants has more than one child
                # If so, it's a convergent synthesis
                child_count = 0
                for child in node.get("children", []):
                    if child["type"] == "mol" and len(child.get("children", [])) > 0:
                        child_count += 1

                if child_count > 1:
                    is_linear = False
                    print("Found convergent step")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # A linear fragment assembly should have multiple fragment couplings and be linear
    if fragment_couplings >= 2 and is_linear:
        print(f"Found linear fragment assembly with {fragment_couplings} couplings")
        return True

    return False
