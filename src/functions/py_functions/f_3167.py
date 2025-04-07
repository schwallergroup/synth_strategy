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
    Detects if the synthesis route includes multiple N-alkylation steps.
    """
    n_alkylation_count = 0

    def dfs_traverse(node):
        nonlocal n_alkylation_count

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check if reactants contain an amine
            reactants = reactants_smiles.split(".")
            has_amine = any(
                Chem.MolFromSmiles(r)
                and Chem.MolFromSmiles(r).HasSubstructMatch(Chem.MolFromSmarts("[#7]"))
                for r in reactants
                if r
            )

            # Check if product contains a new C-N bond
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

            if has_amine and product_mol:
                # Check for N-alkylation pattern
                n_alkyl_pattern = Chem.MolFromSmarts("[#7][#6;!$(C=O)]")
                product_matches = len(product_mol.GetSubstructMatches(n_alkyl_pattern))
                reactant_matches = sum(
                    len(mol.GetSubstructMatches(n_alkyl_pattern))
                    for mol in reactant_mols
                    if mol
                )

                if product_matches > reactant_matches:
                    print(
                        f"Found N-alkylation step, count now: {n_alkylation_count + 1}"
                    )
                    n_alkylation_count += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return n_alkylation_count >= 2
