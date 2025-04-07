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
    Detects Boc protection/deprotection sequence for amines.
    """
    has_boc_protection = False
    has_boc_deprotection = False

    def dfs_traverse(node):
        nonlocal has_boc_protection, has_boc_deprotection

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Boc group pattern
                boc_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])-[#7]")

                # Check for Boc deprotection: Boc in reactants but not in product
                if any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(boc_pattern)
                    for r in reactants
                    if r
                ):
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and not product_mol.HasSubstructMatch(boc_pattern):
                        has_boc_deprotection = True
                        print("Found Boc deprotection")

                # Check for Boc protection: Boc not in reactants but in product
                if not any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(boc_pattern)
                    for r in reactants
                    if r
                ):
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(boc_pattern):
                        has_boc_protection = True
                        print("Found Boc protection")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return has_boc_protection or has_boc_deprotection
