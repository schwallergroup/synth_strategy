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
    """Check if the route involves convergent synthesis (multiple branches converging)."""
    # Count the number of branches at each depth
    branches_by_depth = {}

    def count_branches(node, depth=0):
        if node["type"] == "reaction":
            # Count reactants (children) for this reaction
            reactant_count = sum(1 for child in node.get("children", []) if child["type"] == "mol")
            if reactant_count > 1:
                branches_by_depth[depth] = branches_by_depth.get(depth, 0) + 1

        for child in node.get("children", []):
            count_branches(child, depth + 1)

    count_branches(route)

    # If we have multiple branches at any depth, it's convergent
    is_convergent = any(count > 0 for count in branches_by_depth.values())
    if is_convergent:
        print(f"Found convergent synthesis with branches: {branches_by_depth}")
    return is_convergent
