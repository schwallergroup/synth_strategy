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
    Detects late-stage formylation (introduction of aldehyde group to aromatic ring).
    Late stage is defined as occurring in the first half of the synthesis (lower depth).
    """
    formylation_found = False
    max_depth = 0

    # First pass to find max depth
    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)
        for child in node.get("children", []):
            find_max_depth(child, current_depth + 1)

    find_max_depth(route)

    # Second pass to check for late-stage formylation
    def dfs_traverse(node, depth=0):
        nonlocal formylation_found

        # Consider it late stage if in first half of synthesis (lower depth numbers)
        if depth <= max_depth / 2 and node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aromatic ring in reactants
            aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
            aromatic_present = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(aromatic_pattern):
                    aromatic_present = True

            # Check for new aldehyde in product
            aldehyde_pattern = Chem.MolFromSmarts("[cX3]!@[CX3H1](=O)")
            product_mol = Chem.MolFromSmiles(product)

            if aromatic_present and product_mol and product_mol.HasSubstructMatch(aldehyde_pattern):
                formylation_found = True
                print(f"Found late-stage formylation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return formylation_found
