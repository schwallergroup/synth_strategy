#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    Detects a linear synthesis strategy with multiple protection steps.
    """
    # Track protection steps and their depths
    protection_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if not product or not all(reactants):
                return

            # Check for protection reactions
            alcohol_pattern = Chem.MolFromSmarts("[#6]-[OH]")
            protected_pattern = Chem.MolFromSmarts(
                "[#6]-[#8]-[!#1;!#6;!#8]"
            )  # O connected to non-H, non-C, non-O

            for reactant in reactants:
                if reactant.HasSubstructMatch(alcohol_pattern):
                    if product and product.HasSubstructMatch(protected_pattern):
                        print(f"Found protection at depth {depth}")
                        protection_depths.append(depth)

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found multiple protection steps at different depths
    # indicating a linear protection strategy
    return len(protection_depths) >= 2 and len(set(protection_depths)) >= 2
