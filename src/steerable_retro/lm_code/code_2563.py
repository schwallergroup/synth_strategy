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
    This function detects if the synthesis includes N-alkylation of a heterocycle
    with a benzyl or similar halide.
    """
    found_n_alkylation = False

    def dfs_traverse(node):
        nonlocal found_n_alkylation

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Patterns for N-alkylation components
            benzyl_halide_pattern = Chem.MolFromSmarts("c[CH2][Br,Cl,I]")
            nh_heterocycle_pattern = Chem.MolFromSmarts("[nH]")
            n_alkylated_pattern = Chem.MolFromSmarts("n[CH2]c")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and all(r for r in reactant_mols):
                has_benzyl_halide = any(
                    r.HasSubstructMatch(benzyl_halide_pattern) for r in reactant_mols if r
                )
                has_nh_heterocycle = any(
                    r.HasSubstructMatch(nh_heterocycle_pattern) for r in reactant_mols if r
                )
                has_n_alkylated_product = (
                    product_mol.HasSubstructMatch(n_alkylated_pattern) if product_mol else False
                )

                if has_benzyl_halide and has_nh_heterocycle and has_n_alkylated_product:
                    found_n_alkylation = True
                    print("Found N-alkylation of heterocycle")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_n_alkylation
