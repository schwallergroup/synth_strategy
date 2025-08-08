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
    Detects if the synthesis route involves a cyclohexyl group
    that remains present throughout the synthesis.
    """
    steps_with_cyclohexyl = 0
    total_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal steps_with_cyclohexyl, total_steps

        if node["type"] == "reaction":
            total_steps += 1
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check for cyclohexyl group
                cyclohexyl_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#6][#6]1")

                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(cyclohexyl_pattern):
                    steps_with_cyclohexyl += 1
                    print(f"Cyclohexyl group detected in product at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    # Return True if cyclohexyl is present in at least 75% of steps
    return total_steps > 0 and steps_with_cyclohexyl / total_steps >= 0.75
