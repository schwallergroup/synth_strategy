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
    This function detects a synthesis that preserves a bicyclic core structure
    throughout the synthesis.
    """
    bicyclic_steps = 0
    total_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal bicyclic_steps, total_steps

        if node["type"] == "reaction":
            total_steps += 1

            # Extract product
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Check for bicyclic core
            bicyclic_pattern = Chem.MolFromSmarts("[#6]1[#6]2[#6][#6]2[#6][#7]1")
            if product_mol.HasSubstructMatch(bicyclic_pattern):
                bicyclic_steps += 1
                print(f"Detected bicyclic core at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if bicyclic core is preserved throughout most of the synthesis
    core_preserved = (bicyclic_steps > 0) and (bicyclic_steps / total_steps >= 0.75)
    print(f"Strategy detection result: {core_preserved}")
    print(f"Steps with bicyclic core: {bicyclic_steps}/{total_steps}")

    return core_preserved
