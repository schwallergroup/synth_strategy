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
    This function detects a linear synthesis strategy where each reaction builds
    directly on the product of the previous reaction (as opposed to convergent synthesis).
    """
    # Track reaction depths and number of reactants
    reaction_info = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                depth = node.get("metadata", {}).get("depth", -1)

                # Count number of non-trivial reactants (excluding small molecules)
                significant_reactants = 0
                for r in reactants:
                    mol = Chem.MolFromSmiles(r)
                    if mol and mol.GetNumHeavyAtoms() > 3:  # Consider only substantial fragments
                        significant_reactants += 1

                reaction_info.append((depth, significant_reactants))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # Sort by depth
    reaction_info.sort(key=lambda x: x[0])

    # Check if all reactions have only one significant reactant
    # (characteristic of linear synthesis)
    is_linear = all(num_reactants == 1 for _, num_reactants in reaction_info)

    if is_linear:
        print("Detected linear synthesis strategy")
    else:
        print("Detected convergent or mixed synthesis strategy")

    return is_linear
