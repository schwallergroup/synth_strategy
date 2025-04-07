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
    This function detects if the synthesis involves a convergent amide coupling strategy
    where two fragments are joined via amide bond formation.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is an amide formation reaction
            if len(reactants_smiles) >= 2:  # At least two reactants
                try:
                    # Look for carboxylic acid and amine patterns
                    acid_pattern = Chem.MolFromSmarts("[C$(C=O)][OH]")
                    amine_pattern = Chem.MolFromSmarts("[N;!$(NC=O);!$(N=*)]")
                    amide_pattern = Chem.MolFromSmarts("[C$(C=O)][N]")

                    # Check reactants for acid and amine
                    has_acid = False
                    has_amine = False

                    for r_smiles in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r_smiles)
                        if r_mol:
                            if r_mol.HasSubstructMatch(acid_pattern):
                                has_acid = True
                            if r_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True

                    # Check product for amide
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    has_amide = product_mol and product_mol.HasSubstructMatch(
                        amide_pattern
                    )

                    if has_acid and has_amine and has_amide:
                        print("Convergent amide coupling detected")
                        result = True
                except Exception as e:
                    print(f"Error in amide coupling detection: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return result
