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
    This function detects if the synthetic route involves formation of a heterocycle
    in the final step of the synthesis.
    """
    heterocycle_formed = False
    final_step = True

    def dfs_traverse(node):
        nonlocal heterocycle_formed, final_step

        if node["type"] == "reaction" and final_step:
            # Get reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if heterocycle is formed
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Define heterocycle patterns
            heterocycle_patterns = [
                Chem.MolFromSmarts("[n]1[n][c][nH][c]1=[O]"),  # triazolone
                Chem.MolFromSmarts("[n]1[n][c][n][c]1"),  # triazole
                Chem.MolFromSmarts("[o,n]1[c][c][c][c]1"),  # furan/pyrrole
                Chem.MolFromSmarts("[n]1[c][c][n][c]1"),  # pyrimidine
            ]

            # Check if product contains heterocycle not present in reactants
            for pattern in heterocycle_patterns:
                if product_mol.HasSubstructMatch(pattern):
                    has_pattern_in_reactants = False
                    for r_mol in reactant_mols:
                        if r_mol and r_mol.HasSubstructMatch(pattern):
                            has_pattern_in_reactants = True
                            break

                    if not has_pattern_in_reactants:
                        print("Heterocycle formation detected in final step")
                        heterocycle_formed = True
                        break

            final_step = False  # Mark that we've processed the final step

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return heterocycle_formed
