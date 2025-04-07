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
    This function detects a synthetic strategy involving a nitro-aldol condensation
    (Henry reaction) for C-C bond formation.
    """
    found_nitro_aldol = False

    def dfs_traverse(node, depth=0):
        nonlocal found_nitro_aldol

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Convert to RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if len(reactant_mols) >= 2 and product_mol:
                    # Patterns for aldehyde, nitro compound, and nitro-olefin
                    aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
                    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                    nitro_olefin_pattern = Chem.MolFromSmarts("C=C[N+](=O)[O-]")

                    # Check if reactants contain aldehyde and nitro compound
                    has_aldehyde = any(
                        mol.HasSubstructMatch(aldehyde_pattern) for mol in reactant_mols
                    )
                    has_nitro = any(mol.HasSubstructMatch(nitro_pattern) for mol in reactant_mols)

                    # Check if product is a nitro-olefin
                    product_is_nitro_olefin = product_mol.HasSubstructMatch(nitro_olefin_pattern)

                    if has_aldehyde and has_nitro and product_is_nitro_olefin:
                        found_nitro_aldol = True
                        print(f"Nitro-aldol condensation detected at depth {depth}")
            except:
                print("Error processing reaction SMILES")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_nitro_aldol
