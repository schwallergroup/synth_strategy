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
    This function detects if the synthesis follows a predominantly linear strategy
    with a convergent step only at the end.
    """
    linear_steps = 0
    convergent_steps = 0
    convergent_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal linear_steps, convergent_steps

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count valid reactants (non-empty strings)
                valid_reactants = [r for r in reactants if r.strip()]

                if len(valid_reactants) > 1:
                    convergent_steps += 1
                    convergent_depths.append(depth)
                    print(
                        f"Found convergent step at depth {depth} with {len(valid_reactants)} reactants"
                    )
                else:
                    linear_steps += 1
                    print(f"Found linear step at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # If there are no convergent steps, it's not a linear synthesis with late convergence
    if convergent_steps == 0:
        print("No convergent steps found")
        return False

    # Check if there's only one convergent step and it's at a deeper level (late in synthesis)
    if convergent_steps == 1 and convergent_depths[0] > 1:
        print(
            f"Synthesis has one convergent step at depth {convergent_depths[0]} with only linear steps before it"
        )
        return True

    # If there's only one convergent step at depth 0 or 1, it's not a late convergent step
    if convergent_steps == 1 and convergent_depths[0] <= 1:
        print(
            f"Convergent step is at the beginning (depth {convergent_depths[0]}), not a linear synthesis with late convergence"
        )
        return False

    # If there are multiple convergent steps, check if all except one are at the beginning
    if convergent_steps > 1:
        # Sort depths in ascending order
        sorted_depths = sorted(convergent_depths)

        # If all convergent steps except one are at the beginning (depths 0 and 1)
        early_convergent = [d for d in sorted_depths if d <= 1]
        late_convergent = [d for d in sorted_depths if d > 1]

        if len(early_convergent) == convergent_steps - 1 and len(late_convergent) == 1:
            print(
                f"Synthesis has multiple early convergent steps and one late convergent step at depth {late_convergent[0]}"
            )
            return True

    print("Synthesis does not follow a linear strategy with convergent step only at the end")
    return False
