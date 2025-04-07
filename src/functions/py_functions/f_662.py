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
    This function detects a strategy involving secondary amine formation via nucleophilic substitution,
    typically between a primary amine and a benzyl halide.
    """
    found_secondary_amine_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_secondary_amine_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and all(reactants):
                # Check for primary amine
                primary_amine_pattern = Chem.MolFromSmarts("[NH2][c,C]")
                # Check for benzyl halide
                benzyl_halide_pattern = Chem.MolFromSmarts("[c][C][Br,Cl,I]")
                # Check for secondary amine in product
                secondary_amine_pattern = Chem.MolFromSmarts("[c,C][NH][c,C]")

                has_primary_amine = any(
                    mol.HasSubstructMatch(primary_amine_pattern) for mol in reactants
                )
                has_benzyl_halide = any(
                    mol.HasSubstructMatch(benzyl_halide_pattern) for mol in reactants
                )
                has_secondary_amine = (
                    product.HasSubstructMatch(secondary_amine_pattern)
                    if product
                    else False
                )

                if has_primary_amine and has_benzyl_halide and has_secondary_amine:
                    found_secondary_amine_formation = True
                    print(f"Found secondary amine formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_secondary_amine_formation
