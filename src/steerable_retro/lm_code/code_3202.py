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
    This function detects a synthetic strategy involving N-alkylation of a heterocyclic nitrogen.
    """
    has_n_alkylation = False

    def is_n_alkylation(reactants, product):
        # Check if one reactant has NH and another has alkyl halide
        has_nh = False
        has_alkyl_halide = False

        nh_pattern = Chem.MolFromSmarts("[nH,NH]")
        alkyl_halide_pattern = Chem.MolFromSmarts("C-[#35,#53,#17]")

        for reactant in reactants:
            reactant_mol = Chem.MolFromSmiles(reactant)
            if reactant_mol:
                if reactant_mol.HasSubstructMatch(nh_pattern):
                    has_nh = True
                if reactant_mol.HasSubstructMatch(alkyl_halide_pattern):
                    has_alkyl_halide = True

        # Check if product has N-alkyl bond
        product_mol = Chem.MolFromSmiles(product)
        n_alkyl_pattern = Chem.MolFromSmarts("[n,N]-C")

        if product_mol and product_mol.HasSubstructMatch(n_alkyl_pattern):
            return has_nh and has_alkyl_halide

        return False

    def dfs_traverse(node):
        nonlocal has_n_alkylation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if rsmi:
                parts = rsmi.split(">")
                if len(parts) >= 3:
                    reactants = parts[0].split(".")
                    product = parts[2]

                    # Check for N-alkylation
                    if is_n_alkylation(reactants, product):
                        print("Found N-alkylation reaction")
                        has_n_alkylation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_n_alkylation
