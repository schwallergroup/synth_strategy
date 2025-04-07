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
    Detects if the synthesis includes a reductive amination sequence
    (aldehyde/ketone + amine → amine).
    """
    reductive_amination_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal reductive_amination_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for reductive amination pattern
            aldehyde_ketone_pattern = Chem.MolFromSmarts("[C;$(C=O)]")
            amine_pattern = Chem.MolFromSmarts("[N;!$(N=*);!$(NC=O)]")

            # Check if reactants contain aldehyde/ketone and amine
            has_carbonyl = False
            has_amine = False

            for r_smiles in reactants_smiles:
                r_mol = Chem.MolFromSmiles(r_smiles)
                if r_mol:
                    if r_mol.HasSubstructMatch(aldehyde_ketone_pattern):
                        has_carbonyl = True
                    if r_mol.HasSubstructMatch(amine_pattern):
                        has_amine = True

            # Check if product has a new C-N bond
            if has_carbonyl and has_amine:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                    reductive_amination_detected = True
                    print(f"Reductive amination detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return reductive_amination_detected
