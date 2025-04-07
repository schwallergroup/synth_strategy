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
    This function detects if the synthesis involves reduction of an ester to an alcohol.
    """
    found_pattern = False

    def dfs_traverse(node):
        nonlocal found_pattern

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Ester pattern
            ester_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])-[#6]")
            alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")

            has_ester = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(ester_pattern):
                    has_ester = True

            product_mol = Chem.MolFromSmiles(product)
            if has_ester and product_mol and product_mol.HasSubstructMatch(alcohol_pattern):
                found_pattern = True
                print("Ester reduction to alcohol detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_pattern
