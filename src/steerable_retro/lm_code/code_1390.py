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
    Detects if the synthesis route is primarily linear with a late-stage convergent step.
    This is characterized by having mostly single-reactant steps with a multi-reactant step at low depth.
    """
    convergent_steps = []
    linear_steps = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")

            # Count number of reactants
            num_reactants = len(reactants)

            if num_reactants > 1:
                convergent_steps.append((depth, num_reactants))
                print(f"Detected convergent step at depth {depth} with {num_reactants} reactants")
            else:
                linear_steps.append(depth)
                print(f"Detected linear step at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have a pattern of linear synthesis with late convergence
    if not convergent_steps:
        print("No convergent steps found")
        return False

    # Sort convergent steps by depth
    convergent_steps.sort(key=lambda x: x[0])

    # Check if the earliest convergent step is at a low depth (early in retrosynthesis)
    earliest_convergent_depth = convergent_steps[0][0]
    print(f"Earliest convergent step at depth {earliest_convergent_depth}")

    # Count linear steps at higher depths (later in retrosynthesis)
    linear_steps_after_convergence = sum(1 for d in linear_steps if d > earliest_convergent_depth)
    print(f"Linear steps after convergence: {linear_steps_after_convergence}")

    # Calculate total steps and proportion of linear steps
    total_steps = len(convergent_steps) + len(linear_steps)
    linear_proportion = len(linear_steps) / total_steps if total_steps > 0 else 0
    print(f"Total steps: {total_steps}, Linear proportion: {linear_proportion:.2f}")

    # For a linear with late convergence pattern:
    # 1. We need at least one convergent step at a relatively low depth (≤ 2)
    # 2. We need at least one linear step after the convergent step
    # 3. The overall route should have some linear character
    # 4. The route should have a minimum number of steps
    if (
        earliest_convergent_depth <= 2
        and linear_steps_after_convergence >= 1
        and linear_proportion >= 0.15
        and total_steps >= 3
    ):
        print("Confirmed linear synthesis with late convergence pattern")
        return True

    print("Not a linear synthesis with late convergence pattern")
    return False
