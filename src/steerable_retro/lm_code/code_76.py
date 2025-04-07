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
    This function detects N-alkylation reactions, particularly with alkyl halides.
    """
    n_alkylation_detected = False

    def dfs_traverse(node):
        nonlocal n_alkylation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Patterns for N-alkylation
                nh_pattern = Chem.MolFromSmarts("[nH]")  # Pyrazole NH
                alkyl_halide_pattern = Chem.MolFromSmarts(
                    "[#6][#6]([#6])[I,Br,Cl]"
                )  # Alkyl halide (focusing on isopropyl)
                n_alkylated_pattern = Chem.MolFromSmarts(
                    "[n][#6][#6]([#6])[#6,#1]"
                )  # N-alkylated product

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol and product_mol.HasSubstructMatch(n_alkylated_pattern):
                    if any(r and r.HasSubstructMatch(nh_pattern) for r in reactant_mols):
                        if any(
                            r and r.HasSubstructMatch(alkyl_halide_pattern) for r in reactant_mols
                        ):
                            print("Detected N-alkylation with alkyl halide")
                            n_alkylation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return n_alkylation_detected
