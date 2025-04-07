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
    Detects if the synthesis includes N-methylation of a heterocycle.
    """
    n_methylation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal n_methylation_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for N-methylation pattern
            nh_heterocycle_pattern = Chem.MolFromSmarts("[nH]")
            n_methyl_pattern = Chem.MolFromSmarts("[n][CH3]")
            methyl_donor_pattern = Chem.MolFromSmarts("[CH3][I,Br,Cl,O]")

            # Check reactants and products
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)]
            product_mol = Chem.MolFromSmiles(product)

            if product_mol and product_mol.HasSubstructMatch(n_methyl_pattern):
                # Check if any reactant has NH heterocycle
                if any(
                    mol and mol.HasSubstructMatch(nh_heterocycle_pattern) for mol in reactant_mols
                ):
                    # Check if any reactant is a methyl donor
                    if any(
                        mol and mol.HasSubstructMatch(methyl_donor_pattern) for mol in reactant_mols
                    ):
                        print(f"Found N-methylation of heterocycle at depth {depth}")
                        n_methylation_found = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return n_methylation_found
