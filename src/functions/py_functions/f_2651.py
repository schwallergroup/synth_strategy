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
    Detects a strategy involving multiple esterification steps for protecting carboxylic acids.
    """
    esterification_count = 0

    def dfs_traverse(node):
        nonlocal esterification_count

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for esterification pattern
            reactant_mols = [
                Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
            ]
            product_mol = Chem.MolFromSmiles(product)

            if product_mol and any(
                r.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[OH]"))
                for r in reactant_mols
                if r
            ):
                if product_mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[O][C]")):
                    esterification_count += 1
                    print(
                        f"Found esterification reaction, count: {esterification_count}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found multiple esterification steps
    return esterification_count >= 2
