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
    This function detects a strategy involving heterocycle (isoxazole) formation
    in the synthetic route.
    """
    has_heterocycle_formation = False

    # Define SMARTS pattern for isoxazole
    isoxazole_pattern = Chem.MolFromSmarts("c1conc1")

    def dfs_traverse(node, depth=0):
        nonlocal has_heterocycle_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                # Check if isoxazole is present in product but not in reactants
                if product_mol and product_mol.HasSubstructMatch(isoxazole_pattern):
                    isoxazole_in_reactants = any(
                        r and r.HasSubstructMatch(isoxazole_pattern)
                        for r in reactant_mols
                    )

                    if not isoxazole_in_reactants:
                        has_heterocycle_formation = True
                        print(f"Detected isoxazole formation at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_heterocycle_formation
