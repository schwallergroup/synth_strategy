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
    Detects if the synthesis involves O-alkylation of a phenol.
    """
    phenol_alkylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal phenol_alkylation_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                # Check for phenol in reactants
                phenol_pattern = Chem.MolFromSmarts("c[OH]")
                has_phenol = any(
                    mol.HasSubstructMatch(phenol_pattern)
                    for mol in reactant_mols
                    if mol
                )

                # Check for alkyl halide in reactants
                alkyl_halide_pattern = Chem.MolFromSmarts("[#6][Br,I,Cl]")
                has_alkyl_halide = any(
                    mol.HasSubstructMatch(alkyl_halide_pattern)
                    for mol in reactant_mols
                    if mol
                )

                # Check for aryl ether in product
                aryl_ether_pattern = Chem.MolFromSmarts("c[O][#6]")
                has_aryl_ether = product_mol and product_mol.HasSubstructMatch(
                    aryl_ether_pattern
                )

                if (
                    has_phenol
                    and (has_alkyl_halide or len(reactants) > 1)
                    and has_aryl_ether
                ):
                    print(f"Detected phenol O-alkylation at depth {depth}")
                    phenol_alkylation_detected = True
            except:
                print("Error processing reaction SMILES")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return phenol_alkylation_detected
