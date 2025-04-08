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
    Detects if the synthesis route involves tosylation of an alcohol
    as a leaving group installation strategy.
    """
    found_tosylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_tosylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for alcohol in reactants
                alcohol_pattern = Chem.MolFromSmarts("[OH]")
                # Check for tosylate in product
                tosylate_pattern = Chem.MolFromSmarts("[O][S](=O)(=O)[c]")

                has_alcohol = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(alcohol_pattern)
                    for r in reactants
                    if r
                )

                prod_mol = Chem.MolFromSmiles(product) if product else None
                has_tosylate = prod_mol is not None and prod_mol.HasSubstructMatch(tosylate_pattern)

                if has_alcohol and has_tosylate:
                    found_tosylation = True
                    print(f"Found tosylation at depth {depth}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return found_tosylation
