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
    This function detects if the synthetic route involves nitro reduction to amine.
    """
    nitro_reduction_count = 0

    def dfs_traverse(node):
        nonlocal nitro_reduction_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is a nitro reduction reaction
            nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            product_mol = Chem.MolFromSmiles(product_smiles)

            # Check if product has an amine group
            if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                for reactant_smiles in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    # Check if any reactant has a nitro group
                    if reactant_mol and reactant_mol.HasSubstructMatch(nitro_pattern):
                        # Check if the nitro position in reactant corresponds to amine in product
                        # This is a simplified check and might need refinement
                        nitro_reduction_count += 1
                        print(
                            f"Detected nitro reduction: {reactant_smiles} -> {product_smiles}"
                        )
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Total nitro reductions detected: {nitro_reduction_count}")
    return nitro_reduction_count >= 1
