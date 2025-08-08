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
    This function detects a strategy involving late-stage halogenation
    (introduction of a halogen atom in the final steps).
    """
    # Initialize tracking variable
    has_late_halogenation = False

    # SMARTS pattern for halomethyl group
    halomethyl_pattern = "[#6]-[#6]-[#35,#17,#9,#53]"

    def dfs_traverse(node, depth=0):
        nonlocal has_late_halogenation

        if node["type"] == "reaction" and depth <= 1:  # Late stage (low depth)
            # Extract product
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Check for halomethyl group in product
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts(halomethyl_pattern)
            ):
                has_late_halogenation = True
                print(f"Detected late-stage halogenation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if has_late_halogenation:
        print("Detected late-stage halogenation strategy")
    else:
        print("Late-stage halogenation strategy not detected")

    return has_late_halogenation
