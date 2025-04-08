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
    This function detects if the route contains a nitrile reduction step.
    """
    has_nitrile_reduction = False

    def dfs_traverse(node):
        nonlocal has_nitrile_reduction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles)

            if product and all(r for r in reactants):
                # Check for nitrile reduction
                nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                amine_pattern = Chem.MolFromSmarts("[C][N]")

                reactants_have_nitrile = any(
                    r.HasSubstructMatch(nitrile_pattern) for r in reactants if r
                )
                product_has_amine = product.HasSubstructMatch(amine_pattern)
                product_has_nitrile = product.HasSubstructMatch(nitrile_pattern)

                if reactants_have_nitrile and product_has_amine and not product_has_nitrile:
                    has_nitrile_reduction = True
                    print(f"Nitrile reduction detected in reaction: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Nitrile reduction strategy: {has_nitrile_reduction}")
    return has_nitrile_reduction
