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
    This function detects if the synthetic route maintains an aromatic core with
    an ester functionality throughout the synthesis while modifying side chains.
    """
    # Track if aromatic ester is preserved throughout
    aromatic_ester_preserved = True
    steps_with_aromatic_ester = 0
    total_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal aromatic_ester_preserved, steps_with_aromatic_ester, total_steps

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                total_steps += 1
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Create pattern for aromatic ester
                aromatic_ester_pattern = Chem.MolFromSmarts("[c][C](=[O])[O][C]")

                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(aromatic_ester_pattern):
                    steps_with_aromatic_ester += 1
                    print(f"Aromatic ester found in product at depth {depth}")
                else:
                    aromatic_ester_preserved = False
                    print(f"Aromatic ester not found in product at depth {depth}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if aromatic ester is preserved in all steps and we have at least 2 steps
    return (
        aromatic_ester_preserved and total_steps >= 2 and steps_with_aromatic_ester == total_steps
    )
