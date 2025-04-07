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
    This function detects if the synthesis includes a methoxy deprotection step.
    """
    methoxy_deprotection_detected = False

    def dfs_traverse(node):
        nonlocal methoxy_deprotection_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactants contain methoxy on aromatic but product has hydroxyl instead
            methoxy_pattern = Chem.MolFromSmarts("[#6][#8][#6]:[#6]")
            hydroxy_pattern = Chem.MolFromSmarts("[#8][#6]:[#6]")

            has_methoxy_in_reactants = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(methoxy_pattern):
                    has_methoxy_in_reactants = True
                    break

            product_mol = Chem.MolFromSmiles(product)
            has_hydroxy_in_product = product_mol and product_mol.HasSubstructMatch(
                hydroxy_pattern
            )
            has_methoxy_in_product = product_mol and product_mol.HasSubstructMatch(
                methoxy_pattern
            )

            if (
                has_methoxy_in_reactants
                and has_hydroxy_in_product
                and not has_methoxy_in_product
            ):
                methoxy_deprotection_detected = True
                print("Detected methoxy deprotection")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return methoxy_deprotection_detected
