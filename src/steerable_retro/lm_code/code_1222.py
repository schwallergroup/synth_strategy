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
    Detects a strategy involving sequential reductions: first alkene reduction,
    then nitro group reduction.
    """
    # Initialize flags and depths
    has_alkene_reduction = False
    has_nitro_reduction = False
    alkene_reduction_depth = -1
    nitro_reduction_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal has_alkene_reduction, has_nitro_reduction
        nonlocal alkene_reduction_depth, nitro_reduction_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            if not all(reactants) or not product:
                print(f"Warning: Could not parse all molecules in reaction at depth {depth}")
                return

            # Check for alkene reduction
            c_double_bond_pattern = Chem.MolFromSmarts("C=C")
            if any(
                mol.HasSubstructMatch(c_double_bond_pattern) for mol in reactants
            ) and not product.HasSubstructMatch(c_double_bond_pattern):
                has_alkene_reduction = True
                alkene_reduction_depth = depth
                print(f"Detected alkene reduction at depth {depth}")

            # Check for nitro reduction
            nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            has_nitro_reactant = any(mol.HasSubstructMatch(nitro_pattern) for mol in reactants)
            has_amine_product = product.HasSubstructMatch(amine_pattern)

            if (
                has_nitro_reactant
                and has_amine_product
                and not product.HasSubstructMatch(nitro_pattern)
            ):
                has_nitro_reduction = True
                nitro_reduction_depth = depth
                print(f"Detected nitro reduction at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if both reductions are present and in the correct sequence
    both_reductions = has_alkene_reduction and has_nitro_reduction
    correct_sequence = alkene_reduction_depth > nitro_reduction_depth

    result = both_reductions and correct_sequence

    if result:
        print(
            "Detected sequential reduction strategy: alkene reduction followed by nitro reduction"
        )
    else:
        if both_reductions:
            print("Both reductions present but not in the expected sequence")
        else:
            print("Not all required reductions are present")

    return result
