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
    Detects if the synthesis route preserves a cyclopropylmethoxy group throughout the synthesis.
    """
    cyclopropylmethoxy_pattern = Chem.MolFromSmarts("[#6]1[#6][#6]1[#6][#8]")
    steps_with_cyclopropylmethoxy = 0
    total_steps = 0

    def dfs_traverse(node):
        nonlocal steps_with_cyclopropylmethoxy, total_steps

        if node["type"] == "reaction":
            total_steps += 1
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(cyclopropylmethoxy_pattern):
                    steps_with_cyclopropylmethoxy += 1

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Consider it preserved if it's present in at least 80% of the steps
    if total_steps > 0 and steps_with_cyclopropylmethoxy / total_steps >= 0.8:
        print(
            f"Detected cyclopropylmethoxy preservation in {steps_with_cyclopropylmethoxy}/{total_steps} steps"
        )
        return True
    return False
