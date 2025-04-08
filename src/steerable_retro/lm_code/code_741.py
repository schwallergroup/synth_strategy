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
    This function detects if the synthesis includes a late-stage alcohol oxidation.
    Late stage is defined as occurring in the first half of the synthesis (lower depth).
    """
    max_depth = 0
    has_late_oxidation = False

    def dfs_traverse(node, depth=0):
        nonlocal max_depth, has_late_oxidation

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check for alcohol oxidation
            alcohol_pattern = Chem.MolFromSmarts("[#6][#6;!$(C=O)][OH]")
            carbonyl_pattern = Chem.MolFromSmarts("[#6][C;$(C=O)][#8]")

            reactants = [Chem.MolFromSmiles(smi) for smi in reactants_smiles.split(".")]
            product = Chem.MolFromSmiles(product_smiles)

            if (
                any(mol.HasSubstructMatch(alcohol_pattern) for mol in reactants if mol)
                and product
                and product.HasSubstructMatch(carbonyl_pattern)
            ):
                if depth <= 1:  # Late stage (near the target)
                    has_late_oxidation = True
                    print(f"Late-stage alcohol oxidation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_late_oxidation
