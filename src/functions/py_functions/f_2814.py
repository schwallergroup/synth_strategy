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
    This function detects if the synthesis route contains an N-alkylation
    of a heterocyclic scaffold.
    """
    n_alkylation_found = False

    def dfs_traverse(node):
        nonlocal n_alkylation_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for N-alkylation
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol:
                # Check if reactants contain N-H heterocycle
                nh_heterocycle_pattern = Chem.MolFromSmarts("[nH]")
                n_alkyl_pattern = Chem.MolFromSmarts("[n][C]")

                has_nh = any(
                    mol and mol.HasSubstructMatch(nh_heterocycle_pattern)
                    for mol in reactant_mols
                )
                has_n_alkyl = product_mol.HasSubstructMatch(n_alkyl_pattern)

                if has_nh and has_n_alkyl:
                    print("N-alkylation of heterocycle detected")
                    n_alkylation_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return n_alkylation_found
