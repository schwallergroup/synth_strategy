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
    This function detects if the synthesis follows a convergent approach
    by counting the number of fragments combined in the synthesis.
    """
    max_fragments = 0

    def count_fragments(rsmi):
        reactants = rsmi.split(">")[0].split(".")
        return len(reactants)

    def dfs_traverse(node):
        nonlocal max_fragments

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            fragments = count_fragments(node["metadata"]["rsmi"])
            if fragments > max_fragments:
                max_fragments = fragments

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Maximum fragments in a single reaction: {max_fragments}")
    return max_fragments >= 2  # Consider convergent if at least 2 fragments are combined
