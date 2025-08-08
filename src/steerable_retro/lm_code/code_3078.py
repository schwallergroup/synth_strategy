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
    This function detects a linear synthesis route that incorporates a fluorophenyl group
    throughout the synthesis.
    """
    # Track if we found fluorophenyl in all steps
    steps_with_fluorophenyl = set()
    total_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal total_steps

        if node["type"] == "reaction":
            total_steps += 1

            # Extract product from reaction SMILES
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            product = parts[-1]

            # Check for fluorophenyl group in product
            fluorophenyl_pattern = Chem.MolFromSmarts("c[F]")
            try:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(fluorophenyl_pattern):
                    steps_with_fluorophenyl.add(depth)
                    print(f"Found fluorophenyl group at depth {depth}")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # The strategy is present if fluorophenyl is found in all steps
    strategy_present = len(steps_with_fluorophenyl) == total_steps and total_steps > 0
    print(f"Linear synthesis with fluorophenyl strategy detected: {strategy_present}")
    return strategy_present
