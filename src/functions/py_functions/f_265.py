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
    This function detects a linear synthesis strategy with heterocycle formation.
    """
    linear_synthesis = True
    heterocycle_formation = False

    def dfs_traverse(node):
        nonlocal linear_synthesis, heterocycle_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if more than 2 reactants (would indicate convergent synthesis)
            if len(reactants_smiles) > 2:
                linear_synthesis = False
                print("Convergent synthesis detected (more than 2 reactants)")

            # Check for heterocycle formation
            heterocycle_patterns = [
                Chem.MolFromSmarts("[o;r5]1[c;r5][n;r5][c;r5][c;r5]1"),  # oxazole
                Chem.MolFromSmarts("[n;r5]1[c;r5][c;r5][c;r5][c;r5]1"),  # pyrrole
                Chem.MolFromSmarts(
                    "[n;r6]1[c;r6][c;r6][c;r6][c;r6][c;r6]1"
                ),  # pyridine
            ]

            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol:
                for pattern in heterocycle_patterns:
                    if pattern and product_mol.HasSubstructMatch(pattern):
                        reactants_have_pattern = any(
                            r and r.HasSubstructMatch(pattern)
                            for r in reactants_mols
                            if r
                        )
                        if not reactants_have_pattern:
                            heterocycle_formation = True
                            print("Heterocycle formation detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return linear_synthesis and heterocycle_formation
