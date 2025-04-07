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
    Detects aryl C-N bond formation (likely Buchwald-Hartwig type coupling)
    between an aryl halide and an amine.
    """
    has_aryl_cn_coupling = False

    def dfs_traverse(node):
        nonlocal has_aryl_cn_coupling

        if node.get("type") == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl halide
                aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Br,Cl,I,F]")

                # Check for amine
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                # Check for aryl-N bond in product
                aryl_n_pattern = Chem.MolFromSmarts("[c]-[NH]")

                # Check if reactants contain aryl halide and amine, and product contains aryl-N bond
                has_aryl_halide = any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(aryl_halide_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                )
                has_amine = any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(amine_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                )

                product_mol = Chem.MolFromSmiles(product)
                has_aryl_n = product_mol and product_mol.HasSubstructMatch(
                    aryl_n_pattern
                )

                if has_aryl_halide and has_amine and has_aryl_n:
                    has_aryl_cn_coupling = True
                    print("Found aryl C-N bond formation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_aryl_cn_coupling
