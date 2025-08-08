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
    This function detects if the synthesis follows a linear strategy with a final convergent step.
    """
    # Track synthesis pattern
    is_final_step_convergent = False
    all_other_steps_linear = True

    def dfs_traverse(node):
        nonlocal is_final_step_convergent, all_other_steps_linear

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count non-empty reactants
                reactant_count = sum(1 for r in reactants if r.strip())

                # Check if it's the final step (depth 0)
                if node.get("depth", 0) == 0:
                    if reactant_count >= 2:
                        is_final_step_convergent = True
                        print("Final step is convergent with", reactant_count, "reactants")
                else:
                    # For other steps, check if they're linear
                    if reactant_count > 2:
                        all_other_steps_linear = False
                        print("Found non-linear step at depth", node.get("depth", "unknown"))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return is_final_step_convergent and all_other_steps_linear
