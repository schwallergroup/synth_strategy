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
    Detects if the synthesis involves a Heck coupling reaction (C=C bond formation
    between aryl halide and alkene).
    """
    heck_coupling_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal heck_coupling_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl halide and alkene in reactants, and C=C-aryl in product
            # or the reverse in retrosynthesis
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product)

            if product_mol and all(reactant_mols):
                # For retrosynthesis: product has C=C-aryl, reactants have aryl-X and alkene
                aryl_halide_pattern = Chem.MolFromSmarts("[c][Cl,Br,I]")
                alkene_pattern = Chem.MolFromSmarts("[C]=[C]")
                conjugated_pattern = Chem.MolFromSmarts("[c]/[C]=[C]/[C,c]")

                has_aryl_halide = any(
                    mol.HasSubstructMatch(aryl_halide_pattern) for mol in reactant_mols
                )
                has_alkene = any(
                    mol.HasSubstructMatch(alkene_pattern) for mol in reactant_mols
                )
                has_conjugated = product_mol.HasSubstructMatch(conjugated_pattern)

                if has_aryl_halide and has_alkene and has_conjugated:
                    heck_coupling_detected = True
                    print(f"Heck coupling detected at depth {depth}")

                # For forward synthesis: reactant has C=C-aryl, products have disconnected parts
                if not heck_coupling_detected:
                    has_product_aryl_halide = product_mol.HasSubstructMatch(
                        aryl_halide_pattern
                    )
                    has_product_alkene = product_mol.HasSubstructMatch(alkene_pattern)
                    has_reactant_conjugated = any(
                        mol.HasSubstructMatch(conjugated_pattern)
                        for mol in reactant_mols
                    )

                    if (
                        has_product_aryl_halide
                        and has_product_alkene
                        and has_reactant_conjugated
                    ):
                        heck_coupling_detected = True
                        print(f"Heck coupling (reverse) detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return heck_coupling_detected
