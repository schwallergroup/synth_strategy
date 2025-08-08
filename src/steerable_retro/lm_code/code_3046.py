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
    This function detects a linear fragment assembly strategy where fragments
    are added sequentially rather than in a convergent manner.
    """
    # Track the number of reactants at each step
    reactant_counts = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count the number of reactants
                reactant_counts.append(len(reactants))

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if most reactions have 1-2 reactants (linear strategy)
    # and at least one key step has 2 reactants (fragment coupling)
    if reactant_counts:
        linear_steps = sum(1 for count in reactant_counts if count <= 2)
        has_coupling = any(count == 2 for count in reactant_counts)

        is_linear = (linear_steps / len(reactant_counts) >= 0.7) and has_coupling

        if is_linear:
            print(
                f"Found linear fragment assembly strategy with {linear_steps}/{len(reactant_counts)} linear steps"
            )

        return is_linear

    return False
