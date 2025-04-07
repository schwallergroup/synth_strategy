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
    Detects if the synthesis involves O-alkylation of a phenol with an alkyl halide
    to form an ether linkage.
    """
    has_o_alkylation = False

    def dfs_traverse(node):
        nonlocal has_o_alkylation

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Look for phenol and alkyl halide patterns
            phenol_pattern = Chem.MolFromSmarts("c[OH]")
            alkyl_halide_pattern = Chem.MolFromSmarts("[CX4][Br,Cl,I,F]")

            phenol_found = False
            alkyl_halide_found = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(phenol_pattern):
                        phenol_found = True
                    if mol.HasSubstructMatch(alkyl_halide_pattern):
                        alkyl_halide_found = True

            # Check if product has new C-O-C linkage
            if phenol_found and alkyl_halide_found:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    ether_pattern = Chem.MolFromSmarts("cOC")
                    if product_mol.HasSubstructMatch(ether_pattern):
                        print("Found O-alkylation to form ether")
                        has_o_alkylation = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_o_alkylation
