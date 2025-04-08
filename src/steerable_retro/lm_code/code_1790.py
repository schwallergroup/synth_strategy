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
    This function detects if the synthetic route involves N-alkylation of a
    spirocyclic lactam using a benzyl halide.
    """
    spirocyclic_alkylation_found = False

    def dfs_traverse(node):
        nonlocal spirocyclic_alkylation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for spirocyclic lactam
            spirocyclic_lactam_pattern = Chem.MolFromSmarts(
                "[N]1[C](=[O])[C]2([CH2][CH2]1)[CH2][CH2][CH2][CH2][CH2]2"
            )

            # Check for benzyl halide
            benzyl_halide_pattern = Chem.MolFromSmarts("[c][CH2][Br,Cl,I,F]")

            # Check for N-alkylated product
            n_alkylated_pattern = Chem.MolFromSmarts("[N]([C](=[O]))[CH2][c]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and product_mol.HasSubstructMatch(n_alkylated_pattern):
                # Check if reactants include spirocyclic lactam and benzyl halide
                has_spirocyclic = False
                has_benzyl_halide = False

                for r_mol in reactant_mols:
                    if r_mol:
                        if r_mol.HasSubstructMatch(spirocyclic_lactam_pattern):
                            has_spirocyclic = True
                        if r_mol.HasSubstructMatch(benzyl_halide_pattern):
                            has_benzyl_halide = True

                if has_spirocyclic and has_benzyl_halide:
                    spirocyclic_alkylation_found = True
                    print(f"Found spirocyclic lactam alkylation: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    print(f"Spirocyclic lactam alkylation strategy detected: {spirocyclic_alkylation_found}")
    return spirocyclic_alkylation_found
