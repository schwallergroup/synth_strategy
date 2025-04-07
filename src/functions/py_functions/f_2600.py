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
    This function detects if the synthesis involves early benzylic functionalization
    (bromination or other modification of a benzylic position in the first half of the synthesis).
    """
    has_benzylic_functionalization = False
    max_depth = 0

    # First pass to determine the maximum depth
    def get_max_depth(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        for child in node.get("children", []):
            get_max_depth(child, depth + 1)

    # Second pass to check for benzylic functionalization
    def dfs_traverse(node, depth=0):
        nonlocal has_benzylic_functionalization

        if (
            node["type"] == "reaction" and depth >= max_depth / 2
        ):  # Early stage (second half of max_depth)
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product_mol is not None:
                # Check for benzylic functionalization patterns
                benzylic_br_pattern = Chem.MolFromSmarts("c[C][Br]")  # Benzylic bromide
                benzylic_n_pattern = Chem.MolFromSmarts("c[C][N]")  # Benzylic amine
                benzylic_o_pattern = Chem.MolFromSmarts(
                    "c[C][O]"
                )  # Benzylic alcohol/ether

                if (
                    product_mol.HasSubstructMatch(benzylic_br_pattern)
                    or product_mol.HasSubstructMatch(benzylic_n_pattern)
                    or product_mol.HasSubstructMatch(benzylic_o_pattern)
                ):
                    print(f"Benzylic functionalization detected at depth {depth}")
                    has_benzylic_functionalization = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Get the maximum depth first
    get_max_depth(route)

    # Start traversal from the root
    dfs_traverse(route)
    return has_benzylic_functionalization
