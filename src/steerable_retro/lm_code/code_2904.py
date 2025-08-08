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
    Detects a linear fragment assembly strategy where multiple fragments
    are sequentially added without convergent steps.
    """
    # Initialize tracking variables
    fragment_count = 0
    is_linear = True
    reaction_sequence = []

    def dfs_traverse(node):
        nonlocal fragment_count, is_linear

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                # Count number of reactants
                if len(reactants) > 1:
                    fragment_count += len(reactants) - 1
                    reaction_sequence.append(len(reactants))
                    print(f"Reaction with {len(reactants)} reactants: {rsmi}")

                # Check if any reaction has more than 2 reactants (might indicate convergent synthesis)
                if len(reactants) > 2:
                    is_linear = False
                    print(f"Detected potentially convergent step with {len(reactants)} reactants")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present (at least 3 fragments added in a linear fashion)
    strategy_present = fragment_count >= 2 and is_linear

    if strategy_present:
        print(f"Detected linear fragment assembly with {fragment_count+1} total fragments")
    else:
        print(
            f"Linear fragment assembly not detected or not significant (fragments: {fragment_count+1}, linear: {is_linear})"
        )

    return strategy_present
