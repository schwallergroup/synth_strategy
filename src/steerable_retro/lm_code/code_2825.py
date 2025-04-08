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
    Detects synthesis routes that include N-alkylation with a complex aromatic fragment.
    """
    found_n_alkylation = False

    def dfs_traverse(node):
        nonlocal found_n_alkylation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for N-alkylation
            nh_pattern = Chem.MolFromSmarts("[nH]")
            n_alkyl_pattern = Chem.MolFromSmarts("n[CH2]c")
            aromatic_pattern = Chem.MolFromSmarts("c1ccc2ccccc2c1")  # Naphthalene or similar

            product_mol = Chem.MolFromSmiles(product_smiles)
            reactant_mols = [
                Chem.MolFromSmiles(r) for r in reactants_smiles if Chem.MolFromSmiles(r)
            ]

            if product_mol and product_mol.HasSubstructMatch(n_alkyl_pattern):
                # Check if any reactant has NH and another has aromatic fragment
                has_nh = any(
                    r and r.HasSubstructMatch(nh_pattern) for r in reactant_mols if r is not None
                )
                has_aromatic = any(
                    r and r.HasSubstructMatch(aromatic_pattern)
                    for r in reactant_mols
                    if r is not None
                )

                if has_nh and has_aromatic:
                    print("Found N-alkylation with complex aromatic fragment")
                    found_n_alkylation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_n_alkylation
