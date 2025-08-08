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
    This function detects a linear synthesis strategy (as opposed to convergent)
    by checking if each reaction has at most 2 reactants.
    """
    # Track if all reactions are linear (at most 2 reactants)
    is_linear = True
    # Count total reactions
    reaction_count = 0

    def dfs_traverse(node):
        nonlocal is_linear, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1

            # Get reactants
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if rsmi:
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Count non-empty reactants
                reactant_count = sum(1 for r in reactants_smiles if r)

                # If more than 2 reactants, it's not a linear synthesis
                if reactant_count > 2:
                    is_linear = False
                    print(f"Found non-linear step with {reactant_count} reactants")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we have at least 2 reactions and all are linear
    return reaction_count >= 2 and is_linear
