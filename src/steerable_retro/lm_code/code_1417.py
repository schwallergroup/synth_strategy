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
    This function detects if the synthetic route involves aryl bromide to aldehyde
    transformation as a key step.
    """
    aryl_br_to_aldehyde = False

    def dfs_traverse(node, depth=0):
        nonlocal aryl_br_to_aldehyde

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl bromide in reactants
            aryl_br_pattern = Chem.MolFromSmarts("c[Br]")
            aryl_cho_pattern = Chem.MolFromSmarts("c[CH]=O")

            has_aryl_br = any(
                Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).HasSubstructMatch(aryl_br_pattern)
                for r in reactants
            )
            has_aryl_cho = Chem.MolFromSmiles(product) and Chem.MolFromSmiles(
                product
            ).HasSubstructMatch(aryl_cho_pattern)

            if has_aryl_br and has_aryl_cho:
                aryl_br_to_aldehyde = True
                print(f"Aryl bromide to aldehyde transformation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return aryl_br_to_aldehyde
