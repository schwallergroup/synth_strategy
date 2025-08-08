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
    This function detects a linear fragment assembly strategy by checking if most reactions
    involve combining exactly two fragments.
    """
    reaction_count = 0
    two_fragment_reactions = 0

    def dfs_traverse(node):
        nonlocal reaction_count, two_fragment_reactions

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]

                # Count fragments in reactants
                fragment_count = len(reactants_part.split("."))
                reaction_count += 1

                if fragment_count == 2:
                    two_fragment_reactions += 1
                    print(f"Detected two-fragment reaction: {reactants_part}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If at least 60% of reactions involve exactly two fragments, consider it linear assembly
    return reaction_count > 0 and (two_fragment_reactions / reaction_count) >= 0.6
