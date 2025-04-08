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
    boc_deprotection_detected = False

    def dfs_traverse(node):
        nonlocal boc_deprotection_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactants contain Boc group but product doesn't
            boc_pattern = Chem.MolFromSmarts("[#6]([#6])([#6])([#6])[#8][#6](=[#8])[#7]")

            has_boc_in_reactants = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(boc_pattern):
                    has_boc_in_reactants = True
                    break

            product_mol = Chem.MolFromSmiles(product)
            has_boc_in_product = product_mol and product_mol.HasSubstructMatch(boc_pattern)

            if has_boc_in_reactants and not has_boc_in_product:
                boc_deprotection_detected = True
                print("Detected Boc deprotection")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return boc_deprotection_detected
