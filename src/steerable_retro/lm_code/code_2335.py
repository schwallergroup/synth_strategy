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
    Detects if the synthesis uses amide bond formation to connect fragments.
    """
    amide_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for carboxylic acid in reactants
            carboxylic_acid_pattern = Chem.MolFromSmarts("[#6](=[#8])[#8;H1]")

            # Check for amine in reactants
            amine_pattern = Chem.MolFromSmarts("[#7;H2]")

            # Check for amide in product
            amide_pattern = Chem.MolFromSmarts("[#6](=[#8])[#7]")

            has_acid = False
            has_amine = False
            for r in reactants_smiles:
                mol = Chem.MolFromSmiles(r)
                if mol:
                    if mol.HasSubstructMatch(carboxylic_acid_pattern):
                        has_acid = True
                    if mol.HasSubstructMatch(amine_pattern):
                        has_amine = True

            product_mol = Chem.MolFromSmiles(product_smiles)
            has_amide = product_mol and product_mol.HasSubstructMatch(amide_pattern)

            if has_acid and has_amine and has_amide:
                amide_formation_detected = True
                print(f"Amide formation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Amide formation strategy detected: {amide_formation_detected}")
    return amide_formation_detected
