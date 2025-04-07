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
    This function detects if the synthesis includes a Boc deprotection step.
    """
    has_boc_deprotection = False

    def dfs_traverse(node):
        nonlocal has_boc_deprotection

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc group in reactant
            boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)[N]")

            # Check for free amine in product
            free_amine_pattern = Chem.MolFromSmarts("[NH]")

            has_boc = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(boc_pattern):
                    has_boc = True

            product_mol = Chem.MolFromSmiles(product)
            has_free_amine = product_mol and product_mol.HasSubstructMatch(free_amine_pattern)

            if has_boc and has_free_amine and not product_mol.HasSubstructMatch(boc_pattern):
                has_boc_deprotection = True
                print("Detected Boc deprotection strategy")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_boc_deprotection
