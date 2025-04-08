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
    This function detects if the synthesis uses benzyl groups for phenol protection.
    """
    has_benzyl_protection = False

    def dfs_traverse(node):
        nonlocal has_benzyl_protection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for benzyl bromide or similar in reactants
            benzyl_halide_pattern = Chem.MolFromSmarts("c1ccccc1C[Br,Cl,I]")
            phenol_pattern = Chem.MolFromSmarts("[c][OH]")
            benzyl_ether_pattern = Chem.MolFromSmarts("c1ccccc1C[O][c]")

            reactant_has_benzyl_halide = False
            reactant_has_phenol = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(benzyl_halide_pattern):
                            reactant_has_benzyl_halide = True
                        if mol.HasSubstructMatch(phenol_pattern):
                            reactant_has_phenol = True
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product)
                product_has_benzyl_ether = product_mol and product_mol.HasSubstructMatch(
                    benzyl_ether_pattern
                )
            except:
                product_has_benzyl_ether = False

            if reactant_has_benzyl_halide and reactant_has_phenol and product_has_benzyl_ether:
                has_benzyl_protection = True
                print(f"Found benzyl protection: {rsmi}")

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Uses benzyl protection strategy: {has_benzyl_protection}")
    return has_benzyl_protection
