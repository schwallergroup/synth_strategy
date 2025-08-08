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
    This function detects if the synthetic route follows a linear build-up strategy
    rather than a convergent approach.
    """
    # Track the number of reactants at each step
    step_reactant_counts = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactant_count = len(reactants_part.split("."))

                # Store depth and reactant count
                step_reactant_counts.append((depth, reactant_count))

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Sort by depth (ascending)
    step_reactant_counts.sort(key=lambda x: x[0])

    # Check if most steps have only 1-2 reactants (linear synthesis)
    # and no step has more than 3 reactants
    if step_reactant_counts:
        linear_steps = sum(1 for _, count in step_reactant_counts if count <= 2)
        total_steps = len(step_reactant_counts)
        max_reactants = max(count for _, count in step_reactant_counts)

        is_linear = (linear_steps / total_steps >= 0.75) and (max_reactants <= 3)

        print(f"Linear steps: {linear_steps}/{total_steps}, Max reactants: {max_reactants}")
        print(f"Is linear synthesis: {is_linear}")

        return is_linear

    return False
