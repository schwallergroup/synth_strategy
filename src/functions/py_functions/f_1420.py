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
    This function detects a strategy involving benzyl protection and deprotection of a phenolic OH.
    """
    benzylation_depths = []
    debenzylation_depths = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product)

            # Check for benzylation (phenol + benzyl halide -> benzyl ether)
            phenol_pattern = Chem.MolFromSmarts("[OH]c")
            benzyl_halide_pattern = Chem.MolFromSmarts("[Br,Cl,I][CH2]c1ccccc1")
            benzyl_ether_pattern = Chem.MolFromSmarts("O[CH2]c1ccccc1")

            has_phenol = any(
                mol is not None and mol.HasSubstructMatch(phenol_pattern)
                for mol in reactant_mols
            )
            has_benzyl_halide = any(
                mol is not None and mol.HasSubstructMatch(benzyl_halide_pattern)
                for mol in reactant_mols
            )
            has_benzyl_ether_product = (
                product_mol is not None
                and product_mol.HasSubstructMatch(benzyl_ether_pattern)
            )

            if has_phenol and has_benzyl_halide and has_benzyl_ether_product:
                depth = node.get("depth", 0)
                benzylation_depths.append(depth)
                print(f"Found benzylation at depth {depth}")

            # Check for debenzylation or benzyl exchange
            has_benzyl_ether_reactant = any(
                mol is not None and mol.HasSubstructMatch(benzyl_ether_pattern)
                for mol in reactant_mols
            )

            if has_benzyl_ether_reactant and (has_phenol or has_benzyl_ether_product):
                depth = node.get("depth", 0)
                debenzylation_depths.append(depth)
                print(f"Found debenzylation or benzyl exchange at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have both benzylation and debenzylation/exchange
    if benzylation_depths and debenzylation_depths:
        print("Found benzyl protection-deprotection strategy")
        return True

    return False
