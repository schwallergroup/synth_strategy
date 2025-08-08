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
    This function detects if the synthesis follows a linear strategy (as opposed to convergent).
    """
    step_count = 0
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal step_count, is_linear

        if node["type"] == "reaction":
            step_count += 1

            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # If more than 2 significant reactants, it might be convergent
                significant_reactants = 0
                for r in reactants:
                    # Count only non-trivial reactants (more than 10 atoms)
                    if len(r) > 10:
                        significant_reactants += 1

                if significant_reactants > 2:
                    print(
                        f"Found potential convergent step at depth {depth} with {significant_reactants} significant reactants"
                    )
                    is_linear = False

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return (
        is_linear and step_count >= 5
    )  # At least 5 steps to be considered a significant linear strategy
