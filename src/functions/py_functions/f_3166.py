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
    Detects if the synthesis route includes an epoxide opening step.
    """
    epoxide_opening_found = False

    def dfs_traverse(node):
        nonlocal epoxide_opening_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check if reactants contain an epoxide
            reactants = reactants_smiles.split(".")
            epoxide_pattern = Chem.MolFromSmarts("[#6]1[#8][#6]1")
            has_epoxide = any(
                Chem.MolFromSmiles(r)
                and Chem.MolFromSmiles(r).HasSubstructMatch(epoxide_pattern)
                for r in reactants
                if r
            )

            # Check if product contains alcohol and no epoxide
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None
            has_alcohol = product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts("[#8;H1]")
            )
            no_epoxide = product_mol and not product_mol.HasSubstructMatch(
                epoxide_pattern
            )

            if has_epoxide and has_alcohol and no_epoxide:
                print("Found epoxide opening step")
                epoxide_opening_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return epoxide_opening_found
