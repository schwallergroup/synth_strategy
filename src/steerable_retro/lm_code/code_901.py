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
    This function detects if the route uses a convergent synthesis strategy
    with 3 or more distinct fragments being combined.
    """
    fragment_count = 0
    max_reactants_in_single_reaction = 0

    def dfs_traverse(node):
        nonlocal fragment_count, max_reactants_in_single_reaction

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count number of reactants in this reaction
                num_reactants = len(reactants)
                max_reactants_in_single_reaction = max(
                    max_reactants_in_single_reaction, num_reactants
                )

                # If this reaction combines multiple fragments, increment fragment count
                if num_reactants >= 2:
                    fragment_count += num_reactants - 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # A convergent synthesis with 3+ fragments would have fragment_count >= 2
    # or at least one reaction with 3+ reactants
    result = fragment_count >= 2 or max_reactants_in_single_reaction >= 3
    print(f"Convergent synthesis with 3+ fragments: {result}")
    print(f"Fragment count: {fragment_count}")
    print(f"Max reactants in a single reaction: {max_reactants_in_single_reaction}")

    return result
