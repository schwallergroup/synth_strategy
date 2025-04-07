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
    Detects synthesis routes that form a triazole ring from a nitrile precursor.
    """
    triazole_formed = False
    nitrile_used = False

    def dfs_traverse(node):
        nonlocal triazole_formed, nitrile_used

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check if reactants contain nitrile
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                for r in reactants:
                    if r and r.HasSubstructMatch(nitrile_pattern):
                        nitrile_used = True

                # Check if product contains triazole
                product = Chem.MolFromSmiles(product_smiles)
                triazole_pattern = Chem.MolFromSmarts("[n]1[n][n][c][n]1")
                if product and product.HasSubstructMatch(triazole_pattern):
                    triazole_formed = True

                # Check if this specific reaction forms triazole from nitrile
                if nitrile_used and triazole_formed:
                    print("Detected triazole formation from nitrile")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return triazole_formed and nitrile_used
