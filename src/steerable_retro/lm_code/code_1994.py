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
    This function detects if the synthesis route employs an amine protection-deprotection strategy,
    specifically looking for acetyl protection and subsequent deprotection.
    """
    protection_found = False
    deprotection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal protection_found, deprotection_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for acetyl protection: primary amine to acetamide
            primary_amine_pattern = Chem.MolFromSmarts("[NH2]-[#6]")
            acetamide_pattern = Chem.MolFromSmarts("[NH]-[C](=[O])-[CH3]")

            # Check for protection
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(acetamide_pattern):
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(primary_amine_pattern):
                        print(f"Amine protection detected at depth {depth}")
                        protection_found = True

            # Check for deprotection
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(primary_amine_pattern):
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(acetamide_pattern):
                        print(f"Amine deprotection detected at depth {depth}")
                        deprotection_found = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return protection_found and deprotection_found
