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
    Detects if the synthesis follows a linear strategy (as opposed to convergent).
    """
    # In a linear synthesis, each reaction typically has only one non-commercial reactant
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Count non-commercial reactants (children of type "mol" that aren't in_stock)
            non_commercial_reactants = sum(
                1
                for child in node.get("children", [])
                if child["type"] == "mol" and not child.get("in_stock", False)
            )

            # If more than one non-commercial reactant, it's likely a convergent step
            if non_commercial_reactants > 1:
                print(
                    f"Found convergent step with {non_commercial_reactants} non-commercial reactants"
                )
                is_linear = False

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return is_linear
