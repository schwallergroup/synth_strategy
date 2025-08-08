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
    Detects if the synthesis uses a linear fragment assembly approach rather than convergent.
    """
    reaction_count = 0
    max_reactants_per_step = 0

    def dfs_traverse(node):
        nonlocal reaction_count, max_reactants_per_step

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count number of reactants in this step
                num_reactants = len([r for r in reactants if r.strip()])
                max_reactants_per_step = max(max_reactants_per_step, num_reactants)
                reaction_count += 1
                print(f"Reaction {reaction_count} has {num_reactants} reactants")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Linear synthesis typically has 2 or fewer reactants per step
    # and doesn't converge multiple complex fragments late in the synthesis
    return reaction_count >= 3 and max_reactants_per_step <= 2
