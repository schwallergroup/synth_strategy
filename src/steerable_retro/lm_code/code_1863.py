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
    This function detects a strategy involving N-alkylation of a heterocycle.
    """
    found_n_alkylation = False

    def dfs_traverse(node):
        nonlocal found_n_alkylation

        if node["type"] == "reaction":
            # Extract reaction information
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for N-alkylation
            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # NH-heterocycle pattern (focusing on benzimidazole)
            nh_heterocycle_pattern = Chem.MolFromSmarts("[nH]1cnc2ccccc12")
            # N-alkylated heterocycle pattern
            n_alkylated_pattern = Chem.MolFromSmarts("[n]1([#6][#6])cnc2ccccc12")

            # Check if reactants have NH-heterocycle and product has N-alkylated heterocycle
            reactants_have_nh = any(
                r is not None and r.HasSubstructMatch(nh_heterocycle_pattern)
                for r in reactants_mols
            )
            product_has_n_alkylated = product_mol is not None and product_mol.HasSubstructMatch(
                n_alkylated_pattern
            )

            if reactants_have_nh and product_has_n_alkylated:
                found_n_alkylation = True
                print("Found N-alkylation of heterocycle")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_n_alkylation
