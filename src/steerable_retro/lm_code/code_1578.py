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
    This function detects linear synthesis strategies that maintain stable functional groups
    like trifluoromethyl and heterocycles throughout the synthesis.
    """
    steps_count = 0
    trifluoromethyl_preserved = True
    heterocycle_preserved = True

    def dfs_traverse(node):
        nonlocal steps_count, trifluoromethyl_preserved, heterocycle_preserved

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            steps_count += 1

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for trifluoromethyl group preservation
            if "[C:]([F:])([F:])[F:]" in product and not any(
                "[C:]([F:])([F:])[F:]" in r for r in reactants
            ):
                trifluoromethyl_preserved = False
                print(f"Trifluoromethyl group not preserved in step: {rsmi}")

            # Check for heterocycle preservation (oxazole)
            if "c1oc" in product and not any("c1oc" in r for r in reactants):
                heterocycle_preserved = False
                print(f"Heterocycle not preserved in step: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(
        f"Steps count: {steps_count}, Trifluoromethyl preserved: {trifluoromethyl_preserved}, Heterocycle preserved: {heterocycle_preserved}"
    )

    # Return True if it's a linear synthesis (>3 steps) with preserved stable groups
    return steps_count >= 3 and trifluoromethyl_preserved and heterocycle_preserved
