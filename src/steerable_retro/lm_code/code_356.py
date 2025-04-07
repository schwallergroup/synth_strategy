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
    Detects if the synthesis route uses N-alkylation to form secondary amines
    """
    n_alkylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal n_alkylation_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for N-alkylation pattern: primary amine + alkyl halide → secondary amine
            primary_amine_pattern = Chem.MolFromSmarts("[N;H2]")
            alkyl_halide_pattern = Chem.MolFromSmarts("[C][Br,Cl,I]")
            secondary_amine_pattern = Chem.MolFromSmarts("[N;H1][C]")

            # Check reactants
            has_primary_amine = False
            has_alkyl_halide = False

            for r in reactants:
                r_mol = Chem.MolFromSmiles(r)
                if r_mol:
                    if r_mol.HasSubstructMatch(primary_amine_pattern):
                        has_primary_amine = True
                    if r_mol.HasSubstructMatch(alkyl_halide_pattern):
                        has_alkyl_halide = True

            # Check product
            product_mol = Chem.MolFromSmiles(product)
            has_secondary_amine = product_mol and product_mol.HasSubstructMatch(
                secondary_amine_pattern
            )

            if has_primary_amine and has_alkyl_halide and has_secondary_amine:
                print(f"N-alkylation detected at depth {depth}")
                n_alkylation_detected = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return n_alkylation_detected
