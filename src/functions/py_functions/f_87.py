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
    This function detects if the synthesis involves N-alkylation of a heterocycle,
    particularly indole N-methylation.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            try:
                # Define patterns for indole NH and N-methylated indole
                indole_nh_pattern = Chem.MolFromSmarts("[nH]1c2ccccc2cc1")
                n_methyl_indole_pattern = Chem.MolFromSmarts("[n;$(n(C))]1c2ccccc2cc1")

                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if (
                    reactants_mol
                    and product_mol
                    and reactants_mol.HasSubstructMatch(indole_nh_pattern)
                    and product_mol.HasSubstructMatch(n_methyl_indole_pattern)
                ):
                    print("N-alkylation of indole detected")
                    result = True
            except Exception as e:
                print(f"Error in N-alkylation detection: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return result
