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
    This function detects a nitro reduction to amine transformation in the synthesis route.
    """
    nitro_to_amine_found = False

    def dfs_traverse(node):
        nonlocal nitro_to_amine_found

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitro group in reactants
            nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(nitro_pattern):
                        # Check for amine group in product
                        amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                            # Verify it's a reduction by checking atom mappings
                            if re.search(r"\[N\+:[0-9]+\]", reactant) and re.search(
                                r"\[NH2:[0-9]+\]", product_smiles
                            ):
                                print("Nitro reduction to amine detected")
                                nitro_to_amine_found = True
                except:
                    continue

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return nitro_to_amine_found
