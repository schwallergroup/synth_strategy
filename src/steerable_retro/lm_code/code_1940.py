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
    This function detects if the synthetic route involves an early-stage formylation
    (introduction of an aldehyde group on an aromatic ring).
    """
    formylation_found = False
    formylation_depth = -1
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal formylation_found, formylation_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product_smiles)
            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]

            # Check for formylation (introduction of aldehyde on aromatic)
            aldehyde_pattern = Chem.MolFromSmarts("[c][C;H1]=O")

            if product_mol and not product_mol.HasSubstructMatch(aldehyde_pattern):
                for r_mol in reactants_mols:
                    if r_mol and r_mol.HasSubstructMatch(aldehyde_pattern):
                        formylation_found = True
                        formylation_depth = depth
                        print(f"Formylation found at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it early stage if it's in the second half of the synthesis (depth > max_depth/2)
    return formylation_found and formylation_depth > max_depth / 2
