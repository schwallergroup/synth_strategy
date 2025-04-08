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
    This function detects a synthetic strategy involving phthalimide protection/deprotection
    for primary amines.
    """
    has_phthalimide_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_phthalimide_deprotection

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for phthalimide in reactants
            phthalimide_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#6][#6]1[C](=O)[N][C](=O)")
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            has_phthalimide = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(phthalimide_pattern)
                for r in reactants
                if Chem.MolFromSmiles(r)
            )
            product_mol = Chem.MolFromSmiles(product)
            has_amine = product_mol and product_mol.HasSubstructMatch(amine_pattern)

            if has_phthalimide and has_amine:
                print(f"Detected phthalimide deprotection at depth {depth}")
                has_phthalimide_deprotection = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_phthalimide_deprotection
