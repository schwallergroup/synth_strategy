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
    This function detects if the synthesis follows a linear strategy (vs. convergent).
    """
    is_linear = True
    max_reactants_per_step = 0

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, max_reactants_per_step

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count number of reactants
                num_reactants = len(reactants)
                max_reactants_per_step = max(max_reactants_per_step, num_reactants)

                # If any step has more than 2 reactants, it might not be strictly linear
                if num_reactants > 2:
                    is_linear = False
                    print(f"Found step with {num_reactants} reactants at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # If max reactants per step is consistently low (1-2), it's likely a linear synthesis
    print(f"Maximum reactants per step: {max_reactants_per_step}")
    return is_linear and max_reactants_per_step <= 2
