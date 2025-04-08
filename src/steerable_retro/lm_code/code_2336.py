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
    Detects if the synthesis involves transformation of a nitrile group.
    """
    nitrile_transformation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_transformation_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Nitrile pattern
            nitrile_pattern = Chem.MolFromSmarts("[#6]#[#7]")

            reactants_mol = Chem.MolFromSmiles(reactants_smiles)
            product_mol = Chem.MolFromSmiles(product_smiles)

            if reactants_mol and product_mol:
                reactants_have_nitrile = reactants_mol.HasSubstructMatch(nitrile_pattern)
                product_has_nitrile = product_mol.HasSubstructMatch(nitrile_pattern)

                # Check if this reaction transforms a nitrile
                if reactants_have_nitrile and not product_has_nitrile:
                    nitrile_transformation_detected = True
                    print(f"Nitrile transformation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Nitrile transformation strategy detected: {nitrile_transformation_detected}")
    return nitrile_transformation_detected
