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
    This function detects if the synthesis includes a Suzuki coupling reaction
    (aryl halide + boronic acid → biaryl).
    """
    found_suzuki = False

    def dfs_traverse(node):
        nonlocal found_suzuki

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Patterns for Suzuki coupling components
            aryl_halide_pattern = Chem.MolFromSmarts("c[Cl,Br,I]")
            boronic_acid_pattern = Chem.MolFromSmarts("cB(O)O")
            biaryl_pattern = Chem.MolFromSmarts("c:c")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and all(r for r in reactant_mols):
                has_aryl_halide = any(
                    r.HasSubstructMatch(aryl_halide_pattern) for r in reactant_mols if r
                )
                has_boronic_acid = any(
                    r.HasSubstructMatch(boronic_acid_pattern) for r in reactant_mols if r
                )
                has_biaryl_in_product = (
                    product_mol.HasSubstructMatch(biaryl_pattern) if product_mol else False
                )

                if has_aryl_halide and has_boronic_acid and has_biaryl_in_product:
                    found_suzuki = True
                    print("Found Suzuki coupling reaction")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_suzuki
