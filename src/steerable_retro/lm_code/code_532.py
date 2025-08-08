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
    Detects if the synthesis involves elaboration of a piperidone scaffold
    through multiple functionalization steps.
    """
    # Track piperidone presence at different depths
    piperidone_depths = []
    functionalization_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal piperidone_depths, functionalization_steps

        if node["type"] == "reaction":
            # Extract product
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]
            product = Chem.MolFromSmiles(product_smiles)

            if product is None:
                print(f"Warning: Could not parse product SMILES at depth {depth}")
                return

            # Check for piperidone scaffold
            piperidone_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#7][#6](=[O])1")
            if product.HasSubstructMatch(piperidone_pattern):
                piperidone_depths.append(depth)

                # Check if this is a functionalization step
                if "." in rsmi.split(">")[0]:  # More than one reactant
                    functionalization_steps += 1
                    print(f"Found piperidone functionalization at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have piperidone throughout and multiple functionalization steps
    strategy_found = len(piperidone_depths) >= 2 and functionalization_steps >= 2

    if strategy_found:
        print(
            f"Found piperidone scaffold elaboration strategy with {functionalization_steps} functionalization steps"
        )
    else:
        print(f"Did not find complete piperidone scaffold elaboration strategy")
        print(f"Piperidone found at depths: {piperidone_depths}")
        print(f"Functionalization steps: {functionalization_steps}")

    return strategy_found
