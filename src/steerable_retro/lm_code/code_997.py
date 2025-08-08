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
    Detects multiple C-N bond formations in the synthesis route.
    Looks for at least 2 reactions where a C-N bond is formed.
    """
    cn_bond_formations = 0

    def dfs_traverse(node):
        nonlocal cn_bond_formations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for C-N bond formation
                # This is a simplified check - in practice, you'd need more sophisticated analysis
                if any("N" in r for r in reactants) and "N" in product:
                    # Check if an aryl halide is present (for arylation)
                    if any(re.search(r"c.*Br|c.*Cl|c.*I", r) for r in reactants):
                        cn_bond_formations += 1
                        print(f"Found C-N bond formation: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    result = cn_bond_formations >= 2
    print(f"Multiple C-N bond formations: {result} (count: {cn_bond_formations})")
    return result
