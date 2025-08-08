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
    This function detects if the synthesis route employs a late-stage nitro reduction strategy.
    Late-stage means the reaction occurs at depth 0 or 1 in the synthesis tree.
    """
    nitro_reduction_found = False
    max_depth = 1  # Define what we consider "late-stage"

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_found

        if node["type"] == "reaction" and depth <= max_depth:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro group in reactant
            reactant_mol = Chem.MolFromSmiles(reactants[0])
            nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")

            # Check for amine group in product
            product_mol = Chem.MolFromSmiles(product)
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            if reactant_mol and product_mol:
                if reactant_mol.HasSubstructMatch(nitro_pattern) and product_mol.HasSubstructMatch(
                    amine_pattern
                ):
                    print(f"Found nitro reduction at depth {depth}")
                    nitro_reduction_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return nitro_reduction_found
