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
    Detects if the synthesis includes a nitrile → carboxylic acid → amide
    transformation sequence.
    """
    # Track if we've seen each transformation
    nitrile_to_acid_found = False
    acid_to_amide_found = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_to_acid_found, acid_to_amide_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitrile to carboxylic acid conversion
            nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
            acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")

            # Check for carboxylic acid to amide conversion
            amide_pattern = Chem.MolFromSmarts("[N][C](=[O])")
            amine_pattern = Chem.MolFromSmarts("[N;H2]")

            # Check reactants and products
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)]
            product_mol = Chem.MolFromSmiles(product)

            if product_mol:
                # Check for nitrile to acid conversion
                if any(
                    mol and mol.HasSubstructMatch(nitrile_pattern) for mol in reactant_mols
                ) and product_mol.HasSubstructMatch(acid_pattern):
                    print(f"Found nitrile to carboxylic acid conversion at depth {depth}")
                    nitrile_to_acid_found = True

                # Check for acid to amide conversion
                if (
                    any(mol and mol.HasSubstructMatch(acid_pattern) for mol in reactant_mols)
                    and any(mol and mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols)
                    and product_mol.HasSubstructMatch(amide_pattern)
                ):
                    print(f"Found carboxylic acid to amide conversion at depth {depth}")
                    acid_to_amide_found = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both transformations are found
    return nitrile_to_acid_found and acid_to_amide_found
