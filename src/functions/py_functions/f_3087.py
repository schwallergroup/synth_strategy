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
    This function detects a synthetic strategy involving Boc deprotection
    of an amine functional group.
    """
    boc_deprotection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_deprotection_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc deprotection
            # Boc protected amine pattern
            boc_amine_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])-[#7]")
            # Primary amine pattern
            primary_amine_pattern = Chem.MolFromSmarts("[#7H2]")

            reactant_has_boc = Chem.MolFromSmiles(
                product
            ) is not None and Chem.MolFromSmiles(product).HasSubstructMatch(
                boc_amine_pattern
            )

            product_has_amine = any(
                Chem.MolFromSmiles(r) is not None
                and Chem.MolFromSmiles(r).HasSubstructMatch(primary_amine_pattern)
                for r in reactants
                if Chem.MolFromSmiles(r)
            )

            if reactant_has_boc and product_has_amine:
                print(f"Found Boc deprotection at depth {depth}")
                boc_deprotection_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Boc deprotection strategy detected: {boc_deprotection_found}")
    return boc_deprotection_found
