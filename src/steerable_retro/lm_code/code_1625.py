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
    Detects a strategy where a fluorophenyl group is introduced in the late stage of synthesis
    via a coupling reaction.
    """
    fluorophenyl_introduction = None

    def dfs_traverse(node):
        nonlocal fluorophenyl_introduction

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                reaction_smiles = node["metadata"]["rsmi"]
                depth = node.get("depth", 0)

                # Check if this is a late stage reaction (depth close to 0)
                if depth <= 1:  # Considering depth 0 or 1 as late stage
                    reactants = reaction_smiles.split(">")[0].split(".")

                    # Look for fluorophenyl pattern in reactants
                    fluorophenyl_pattern = re.compile(r"[cC]1\([F]\)[cC][cC][cC][cC][cC]1")
                    for reactant in reactants:
                        if "OB(O)" in reactant and "F" in reactant:
                            # Check if it's a boronic acid with fluorine
                            fluorophenyl_introduction = depth
                            print(f"Found late-stage fluorophenyl introduction at depth {depth}")
                            break

        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we found a late-stage fluorophenyl introduction
    if fluorophenyl_introduction is not None:
        print("Confirmed late-stage fluorophenyl introduction strategy")
        return True

    return False
