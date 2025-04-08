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
    This function detects a linear synthesis strategy that incorporates multiple
    distinct heterocyclic motifs (quinolinone, thiazole, piperazine).
    """
    # Initialize tracking variables
    has_quinolinone = False
    has_thiazole = False
    has_piperazine = False
    reaction_count = 0

    # Define SMARTS patterns
    quinolinone_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c]2[c]1[n][c](=[O])[c][c]2")
    thiazole_pattern = Chem.MolFromSmarts("[c]1[s][c][n][c]1")
    piperazine_pattern = Chem.MolFromSmarts("[N]1[C][C][N][C][C]1")

    def dfs_traverse(node):
        nonlocal has_quinolinone, has_thiazole, has_piperazine, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1
            rsmi = node.get("metadata", {}).get("rsmi", "")

            if rsmi:
                product = rsmi.split(">")[-1]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol:
                    # Check for heterocyclic motifs
                    if product_mol.HasSubstructMatch(quinolinone_pattern):
                        print("Detected quinolinone motif")
                        has_quinolinone = True

                    if product_mol.HasSubstructMatch(thiazole_pattern):
                        print("Detected thiazole motif")
                        has_thiazole = True

                    if product_mol.HasSubstructMatch(piperazine_pattern):
                        print("Detected piperazine motif")
                        has_piperazine = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if it's a linear synthesis (at least 4 reactions) with heterocycle diversity
    is_linear = reaction_count >= 4
    has_heterocycle_diversity = has_quinolinone and has_thiazole and has_piperazine

    return is_linear and has_heterocycle_diversity
