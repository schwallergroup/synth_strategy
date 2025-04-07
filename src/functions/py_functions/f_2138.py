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
    Detects a synthesis involving phenol alkylation with a haloalkyl chain.
    """
    has_phenol_alkylation = False

    def dfs_traverse(node):
        nonlocal has_phenol_alkylation

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for phenol and haloalkyl patterns in reactants
            phenol_present = False
            haloalkyl_present = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    phenol_pattern = Chem.MolFromSmarts("[OH][c]")
                    if mol.HasSubstructMatch(phenol_pattern):
                        phenol_present = True

                    haloalkyl_pattern = Chem.MolFromSmarts(
                        "[Br,Cl][CH2][CH2][CH2][Cl,Br]"
                    )
                    if mol.HasSubstructMatch(haloalkyl_pattern):
                        haloalkyl_present = True

            # Check for alkylated phenol in product
            if phenol_present and haloalkyl_present:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    alkylated_phenol_pattern = Chem.MolFromSmarts(
                        "[c][O][CH2][CH2][CH2][Cl,Br]"
                    )
                    if product_mol.HasSubstructMatch(alkylated_phenol_pattern):
                        has_phenol_alkylation = True
                        print(f"Found phenol alkylation in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Phenol alkylation strategy detected: {has_phenol_alkylation}")
    return has_phenol_alkylation
