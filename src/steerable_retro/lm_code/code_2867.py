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
    Detects if the synthetic route follows a linear fragment assembly strategy
    rather than a convergent approach.
    """
    # Track the number of reactants at each step
    reactant_counts = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            # Count non-empty reactants
            count = sum(1 for r in reactants_smiles if r)
            reactant_counts.append(count)

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have reactions
    if not reactant_counts:
        return False

    # Linear assembly typically has 2 reactants per step
    # (one main fragment and one small molecule/reagent)
    linear_steps = sum(1 for count in reactant_counts if count == 2)

    # Strategy is present if most steps have 2 reactants
    strategy_present = linear_steps >= len(reactant_counts) * 0.75
    print(f"Linear fragment assembly strategy detected: {strategy_present}")
    return strategy_present
