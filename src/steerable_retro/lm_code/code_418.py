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
    This function detects if the synthetic route follows a linear synthesis strategy
    without convergent steps (each reaction has only one non-commercial reactant).
    """
    is_linear = True
    reaction_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, reaction_count

        indent = "  " * depth
        if node["type"] == "reaction":
            reaction_count += 1
            non_commercial_reactants = 0

            # Count non-commercial reactants
            for child in node.get("children", []):
                if child["type"] == "mol":
                    # Consider a molecule commercial if it's marked in_stock or is a leaf node
                    is_commercial = (
                        child.get("in_stock", False) or len(child.get("children", [])) == 0
                    )
                    if not is_commercial:
                        non_commercial_reactants += 1
                        print(
                            f"{indent}Found non-commercial reactant: {child.get('smiles', 'No SMILES')}"
                        )

            print(
                f"{indent}Reaction {reaction_count} has {non_commercial_reactants} non-commercial reactants"
            )

            # If more than one non-commercial reactant, it's not linear
            if non_commercial_reactants > 1:
                is_linear = False
                print(
                    f"{indent}Found convergent step with {non_commercial_reactants} non-commercial reactants"
                )

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting traversal to check for linear synthesis strategy")
    dfs_traverse(route)

    print(f"Total reactions: {reaction_count}, Is linear: {is_linear}")
    # Must have at least 2 reactions to be considered a strategy
    return is_linear and reaction_count >= 2
