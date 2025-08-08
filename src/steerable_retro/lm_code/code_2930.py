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
    This function detects a predominantly linear synthesis with late-stage fragment coupling.
    """
    # Track branching factor and fragment couplings by depth
    fragment_couplings = []
    reaction_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal fragment_couplings, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]

                # Check for fragment coupling (multiple reactants)
                if reactants.count(".") > 0:
                    fragment_couplings.append((depth, reactants.count(".") + 1))
                    print(
                        f"Detected fragment coupling at depth {depth} with {reactants.count('.') + 1} fragments"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Analyze the pattern
    if reaction_count > 0:
        # Sort couplings by depth
        fragment_couplings.sort(key=lambda x: x[0])

        # Check if most fragment couplings occur late in the synthesis (low depth)
        late_couplings = [fc for fc in fragment_couplings if fc[0] <= 2]
        early_couplings = [fc for fc in fragment_couplings if fc[0] > 2]

        # If there are more late couplings than early ones, or if the most significant coupling is late
        if len(late_couplings) > len(early_couplings) or (
            fragment_couplings and fragment_couplings[0][0] <= 2 and fragment_couplings[0][1] >= 2
        ):
            return True

    return False
