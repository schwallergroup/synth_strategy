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
    This function detects if the synthesis starts with heterocycle formation
    (specifically benzimidazole) and then proceeds with sequential functionalization.
    """
    benzimidazole_formed = False
    benzimidazole_formation_depth = -1
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal benzimidazole_formed, benzimidazole_formation_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this reaction forms a benzimidazole
            product_mol = Chem.MolFromSmiles(product_smiles)
            benzimidazole_pattern = Chem.MolFromSmarts("c1ccc2[nH]cnc2c1")

            if product_mol and benzimidazole_pattern:
                if product_mol.HasSubstructMatch(benzimidazole_pattern):
                    # Check if reactants don't have benzimidazole
                    has_benzimidazole_in_reactants = False
                    for r_smiles in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r_smiles)
                        if r_mol and r_mol.HasSubstructMatch(benzimidazole_pattern):
                            has_benzimidazole_in_reactants = True
                            break

                    if not has_benzimidazole_in_reactants:
                        benzimidazole_formed = True
                        benzimidazole_formation_depth = depth
                        print(f"Benzimidazole formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if benzimidazole formation is at the beginning (highest depth)
    is_heterocycle_first = (
        benzimidazole_formed and benzimidazole_formation_depth >= max_depth - 1
    )

    print(f"Heterocycle-first strategy detected: {is_heterocycle_first}")
    return is_heterocycle_first
