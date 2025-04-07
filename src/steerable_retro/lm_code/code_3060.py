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
    Detects the use of Boc protection/deprotection strategy for nitrogen-containing groups.
    """
    has_boc_protected = False
    has_deprotection = False

    def dfs_traverse(node):
        nonlocal has_boc_protected, has_deprotection

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for Boc group
                boc_pattern = Chem.MolFromSmarts("[CH3][C]([CH3])([CH3])[O][C](=[O])[N]")
                if mol.HasSubstructMatch(boc_pattern):
                    has_boc_protected = True
                    print("Found Boc-protected intermediate")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc deprotection (reactant has Boc, product doesn't)
            reactant_has_boc = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    boc_pattern = Chem.MolFromSmarts("[CH3][C]([CH3])([CH3])[O][C](=[O])[N]")
                    if mol.HasSubstructMatch(boc_pattern):
                        reactant_has_boc = True
                        break

            product_mol = Chem.MolFromSmiles(product)
            product_has_boc = False
            if product_mol:
                boc_pattern = Chem.MolFromSmarts("[CH3][C]([CH3])([CH3])[O][C](=[O])[N]")
                product_has_boc = product_mol.HasSubstructMatch(boc_pattern)

            if reactant_has_boc and not product_has_boc:
                has_deprotection = True
                print("Found Boc deprotection step")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_boc_protected and has_deprotection
