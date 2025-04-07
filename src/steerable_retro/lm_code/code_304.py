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
    This function detects a synthetic strategy involving N-arylation of a heterocycle.
    """
    found_n_arylation = False

    def dfs_traverse(node):
        nonlocal found_n_arylation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for N-arylation pattern
            if len(reactants) >= 2 and product:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if all(mol is not None for mol in reactant_mols) and product_mol is not None:
                    # Look for aryl halide in reactants
                    aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,Cl,I]")
                    nh_heterocycle_pattern = Chem.MolFromSmarts("[nH]")
                    n_aryl_pattern = Chem.MolFromSmarts("[n][c]")

                    has_aryl_halide = any(
                        len(mol.GetSubstructMatches(aryl_halide_pattern)) > 0
                        for mol in reactant_mols
                    )
                    has_nh_heterocycle = any(
                        len(mol.GetSubstructMatches(nh_heterocycle_pattern)) > 0
                        for mol in reactant_mols
                    )
                    has_n_aryl_product = len(product_mol.GetSubstructMatches(n_aryl_pattern)) > 0

                    if has_aryl_halide and has_nh_heterocycle and has_n_aryl_product:
                        found_n_arylation = True
                        print(f"Detected N-arylation of heterocycle: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"N-arylation of heterocycle detected: {found_n_arylation}")
    return found_n_arylation
