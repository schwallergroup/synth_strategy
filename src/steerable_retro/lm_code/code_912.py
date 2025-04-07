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
    Detects if the synthesis includes the conversion of a nitrile to a lactam through
    reduction and cyclization.
    """
    nitrile_to_lactam_found = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_to_lactam_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitrile in reactants
            nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")

            reactants_have_nitrile = False
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(nitrile_pattern):
                        reactants_have_nitrile = True
                        break
                except:
                    continue

            # Check for lactam in product
            lactam_pattern = Chem.MolFromSmarts("[C]1[C][C][N][C](=[O])[C]1")

            product_has_lactam = False
            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(lactam_pattern):
                    product_has_lactam = True
            except:
                pass

            if reactants_have_nitrile and product_has_lactam:
                nitrile_to_lactam_found = True
                print(f"Nitrile to lactam conversion detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return nitrile_to_lactam_found
