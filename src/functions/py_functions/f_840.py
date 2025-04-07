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
    Detects if the synthesis involves multiple functional group interconversions
    (nitro→amine, ester→acid→acid chloride).
    """
    has_nitro_reduction = False
    has_ester_hydrolysis = False
    has_acid_activation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_nitro_reduction, has_ester_hydrolysis, has_acid_activation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product_smiles)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]

            # Check for nitro reduction
            amine_pattern = Chem.MolFromSmarts("[NX3;H2]")
            nitro_pattern = Chem.MolFromSmarts("[N+](=[O-])[O-]")

            if product_mol.HasSubstructMatch(amine_pattern):
                for r_mol in reactant_mols:
                    if r_mol and r_mol.HasSubstructMatch(nitro_pattern):
                        has_nitro_reduction = True
                        print(f"Found nitro reduction at depth {depth}")
                        break

            # Check for ester hydrolysis
            carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H]")
            ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4]")

            if product_mol.HasSubstructMatch(carboxylic_acid_pattern):
                for r_mol in reactant_mols:
                    if r_mol and r_mol.HasSubstructMatch(ester_pattern):
                        has_ester_hydrolysis = True
                        print(f"Found ester hydrolysis at depth {depth}")
                        break

            # Check for acid activation
            acid_chloride_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[Cl]")

            if product_mol.HasSubstructMatch(acid_chloride_pattern):
                for r_mol in reactant_mols:
                    if r_mol and r_mol.HasSubstructMatch(carboxylic_acid_pattern):
                        has_acid_activation = True
                        print(f"Found acid activation at depth {depth}")
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if the strategy is present
    fg_transformations = sum(
        [has_nitro_reduction, has_ester_hydrolysis, has_acid_activation]
    )
    strategy_present = fg_transformations >= 2

    print(f"Functional group interconversion strategy detected: {strategy_present}")
    print(f"Number of functional group transformations: {fg_transformations}")

    return strategy_present
