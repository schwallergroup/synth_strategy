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
    This function detects if the synthetic route follows a linear strategy without convergent steps.
    """
    max_branching = 0

    def count_reactants(reaction_node):
        if "rsmi" in reaction_node.get("metadata", {}):
            rsmi = reaction_node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            return len(reactants)
        return 0

    def dfs_traverse(node):
        nonlocal max_branching

        if node["type"] == "reaction":
            reactant_count = count_reactants(node)
            max_branching = max(max_branching, reactant_count)
            print(f"Reaction with {reactant_count} reactants found")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # A linear synthesis typically has no more than 2 reactants per step
    is_linear = max_branching <= 2
    print(f"Linear synthesis strategy: {is_linear} (Max reactants per step: {max_branching})")
    return is_linear
