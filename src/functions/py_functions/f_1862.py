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
    This function detects a strategy involving nitro reduction to form amines
    as key intermediates in the synthesis.
    """
    found_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal found_nitro_reduction

        if node["type"] == "reaction":
            # Extract reaction information
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitro reduction
            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Nitro pattern
            nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])-[O-]")
            # Amine pattern
            amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")

            # Check if reactants have nitro groups and product has amines
            reactants_have_nitro = any(
                r is not None and r.HasSubstructMatch(nitro_pattern)
                for r in reactants_mols
            )
            product_has_amine = (
                product_mol is not None and product_mol.HasSubstructMatch(amine_pattern)
            )

            if reactants_have_nitro and product_has_amine:
                found_nitro_reduction = True
                print("Found nitro reduction to amine")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_nitro_reduction
