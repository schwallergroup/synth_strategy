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
    This function detects if the synthesis follows a linear (non-convergent) strategy.
    """
    is_linear = True
    max_fragments_per_reaction = 0

    def dfs_traverse(node):
        nonlocal is_linear, max_fragments_per_reaction

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count number of reactant fragments
                fragment_count = sum(1 for r in reactants if r.strip())
                max_fragments_per_reaction = max(max_fragments_per_reaction, fragment_count)

                # If any reaction has more than 2 fragments, it's potentially convergent
                if fragment_count > 2:
                    is_linear = False
                    print(f"Found potential convergent step with {fragment_count} fragments")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # If max fragments is 2 or less and synthesis is linear
    if is_linear and max_fragments_per_reaction <= 2:
        print("Linear synthesis strategy detected")
        return True
    return False
