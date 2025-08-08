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
    Detects if the synthesis follows a linear pattern where each reaction
    adds or modifies one fragment at a time, rather than convergent synthesis.
    """
    # Track if we've found any convergent steps (more than one complex fragment)
    convergent_step_found = False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_step_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count complex reactants (those with more than 10 atoms)
            complex_reactants = 0
            for r in reactants:
                mol = Chem.MolFromSmiles(r)
                if mol and mol.GetNumAtoms() > 10:
                    complex_reactants += 1

            # If more than one complex reactant, it's a convergent step
            if complex_reactants > 1:
                convergent_step_found = True
                print(
                    f"Found convergent step at depth {depth} with {complex_reactants} complex reactants"
                )

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if it's a linear synthesis (no convergent steps)
    return not convergent_step_found
