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
    Detects a synthetic strategy involving malonate alkylation for C-C bond formation.
    """
    has_malonate_alkylation = False

    def dfs_traverse(node):
        nonlocal has_malonate_alkylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                if reactants and product:
                    # Check for malonate pattern in reactants
                    malonate_pattern = Chem.MolFromSmarts("[CH2][C](=[O])[O][CH2][CH3]")
                    alkyl_halide_pattern = Chem.MolFromSmarts("[Br,Cl,I][CH2]")

                    has_malonate = False
                    has_alkyl_halide = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(malonate_pattern):
                                has_malonate = True
                            if reactant_mol.HasSubstructMatch(alkyl_halide_pattern):
                                has_alkyl_halide = True

                    # Check for alkylated malonate in product
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        alkylated_malonate_pattern = Chem.MolFromSmarts(
                            "[CH]([CH2][*])([C](=[O])[O][CH2][CH3])[C](=[O])[O][CH2][CH3]"
                        )
                        if product_mol.HasSubstructMatch(alkylated_malonate_pattern):
                            if has_malonate and has_alkyl_halide:
                                has_malonate_alkylation = True

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Malonate alkylation strategy: {has_malonate_alkylation}")
    return has_malonate_alkylation
