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
    Detects if the synthesis uses an olefination disconnection strategy
    (breaking α,β-unsaturated carbonyl to ketone/aldehyde).
    """
    olefination_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal olefination_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for olefination disconnection pattern
            unsaturated_carbonyl_pattern = Chem.MolFromSmarts("[C;$(C=CC=O)]")
            carbonyl_pattern = Chem.MolFromSmarts("[C;$(C=O);!$(C=OO)]")
            phosphonate_pattern = Chem.MolFromSmarts("[P;$(P(=O)(O)(O))]")

            # Check if product has α,β-unsaturated carbonyl
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol and product_mol.HasSubstructMatch(unsaturated_carbonyl_pattern):
                # Check if reactants have carbonyl and phosphonate
                has_carbonyl = False
                has_phosphonate = False

                for r_smiles in reactants_smiles:
                    r_mol = Chem.MolFromSmiles(r_smiles)
                    if r_mol:
                        if r_mol.HasSubstructMatch(carbonyl_pattern):
                            has_carbonyl = True
                        if r_mol.HasSubstructMatch(phosphonate_pattern):
                            has_phosphonate = True

                if has_carbonyl and has_phosphonate:
                    olefination_detected = True
                    print(f"Olefination disconnection detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return olefination_detected
