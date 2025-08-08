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
    This function detects if the synthesis follows a linear strategy with sequential
    fragment couplings rather than a convergent approach.
    """
    is_linear = True
    reaction_count = 0

    def dfs_traverse(node):
        nonlocal is_linear, reaction_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            reaction_count += 1
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # If there are more than 2 reactants, it might be a convergent step
            if len(reactants) > 2:
                # Check if the extra reactants are small molecules (potential reagents, not fragments)
                large_fragments = 0
                for reactant in reactants:
                    if len(reactant) > 10:  # Simple heuristic for "large" fragments
                        large_fragments += 1

                if large_fragments > 2:
                    print(
                        f"Potential convergent step detected with {large_fragments} large fragments: {rsmi}"
                    )
                    is_linear = False

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Only consider routes with at least 3 reactions
    if reaction_count >= 3:
        return is_linear
    else:
        return False  # Not enough reactions to determine strategy
