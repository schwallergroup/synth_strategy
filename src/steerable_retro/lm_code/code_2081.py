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
    Detects synthesis routes that convert a nitrile to a heterocycle.
    """
    nitrile_to_heterocycle = False

    def dfs_traverse(node):
        nonlocal nitrile_to_heterocycle

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check for nitrile in reactants
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                nitrile_present = any(r and r.HasSubstructMatch(nitrile_pattern) for r in reactants)

                # Check for heterocycle in product
                product = Chem.MolFromSmiles(product_smiles)
                # Generic heterocycle pattern - this is simplified and might need refinement
                heterocycle_pattern = Chem.MolFromSmarts("*1[n,o,s]*[n,o,s]*1")

                if product and product.HasSubstructMatch(heterocycle_pattern) and nitrile_present:
                    nitrile_to_heterocycle = True
                    print("Detected conversion of nitrile to heterocycle")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitrile_to_heterocycle
