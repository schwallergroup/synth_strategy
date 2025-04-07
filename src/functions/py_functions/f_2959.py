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
    This function detects if the synthetic route contains multiple different types of
    heteroatom bond formations (C-N, N-S, C-O).
    """
    c_n_formations = 0
    n_s_formations = 0
    c_o_formations = 0

    def dfs_traverse(node):
        nonlocal c_n_formations, n_s_formations, c_o_formations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants.split(".")]
                product_mol = Chem.MolFromSmiles(product)

                if not product_mol:
                    return

                # Check for C-N bond formation
                c_n_pattern = Chem.MolFromSmarts("[#6][#7]")
                c_n_count_product = len(product_mol.GetSubstructMatches(c_n_pattern))
                c_n_count_reactants = sum(
                    len(mol.GetSubstructMatches(c_n_pattern))
                    for mol in reactant_mols
                    if mol
                )

                if c_n_count_product > c_n_count_reactants:
                    c_n_formations += 1
                    print(f"Found C-N bond formation: {rsmi}")

                # Check for N-S bond formation
                n_s_pattern = Chem.MolFromSmarts("[#7][#16]")
                n_s_count_product = len(product_mol.GetSubstructMatches(n_s_pattern))
                n_s_count_reactants = sum(
                    len(mol.GetSubstructMatches(n_s_pattern))
                    for mol in reactant_mols
                    if mol
                )

                if n_s_count_product > n_s_count_reactants:
                    n_s_formations += 1
                    print(f"Found N-S bond formation: {rsmi}")

                # Check for C-O bond formation
                c_o_pattern = Chem.MolFromSmarts("[#6][#8]")
                c_o_count_product = len(product_mol.GetSubstructMatches(c_o_pattern))
                c_o_count_reactants = sum(
                    len(mol.GetSubstructMatches(c_o_pattern))
                    for mol in reactant_mols
                    if mol
                )

                if c_o_count_product > c_o_count_reactants:
                    c_o_formations += 1
                    print(f"Found C-O bond formation: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Count different types of heteroatom bond formations
    different_types = (c_n_formations > 0) + (n_s_formations > 0) + (c_o_formations > 0)
    result = different_types >= 2

    print(
        f"Multiple heteroatom bond formations detected: {result} (C-N: {c_n_formations}, N-S: {n_s_formations}, C-O: {c_o_formations})"
    )
    return result
