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
    This function detects if the synthetic route follows a linear strategy without
    convergent steps.
    """
    is_linear = True

    def count_reactants(reaction_node):
        if "rsmi" in reaction_node.get("metadata", {}):
            rsmi = reaction_node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            return len(reactants)
        return 0

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # If a reaction has more than 2 reactants, it might be convergent
            reactant_count = count_reactants(node)
            if reactant_count > 2:
                is_linear = False
                print(f"Found potential convergent step with {reactant_count} reactants")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return is_linear
