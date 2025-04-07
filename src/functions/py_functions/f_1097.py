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
    This function detects functionalization of a triazolopyrimidine core.
    """
    # Track if we've found the core heterocycle and its functionalization
    found_triazolopyrimidine = False
    found_functionalization = False

    # SMARTS for triazolopyrimidine core
    triazolopyrimidine_pattern = Chem.MolFromSmarts(
        "[#6]1:[#7]:[#6]:[#7]:[#6]2:[#7]:[#7]:[#6]:[#6]:12"
    )

    def dfs_traverse(node, depth=0):
        nonlocal found_triazolopyrimidine, found_functionalization

        if node["type"] == "mol":
            if "smiles" in node:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol and mol.HasSubstructMatch(triazolopyrimidine_pattern):
                    found_triazolopyrimidine = True
                    print(f"Found triazolopyrimidine core at depth {depth}")

        elif node["type"] == "reaction":
            if found_triazolopyrimidine and "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                reactants = [
                    Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r
                ]
                product = Chem.MolFromSmiles(product_str) if product_str else None

                # Check if both reactants and product contain the core
                reactants_have_core = any(
                    r and r.HasSubstructMatch(triazolopyrimidine_pattern)
                    for r in reactants
                )
                product_has_core = product and product.HasSubstructMatch(
                    triazolopyrimidine_pattern
                )

                if reactants_have_core and product_has_core:
                    found_functionalization = True
                    print(
                        f"Found functionalization of triazolopyrimidine at depth {depth}"
                    )

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    strategy_present = found_triazolopyrimidine and found_functionalization
    print(f"Heterocycle functionalization strategy detected: {strategy_present}")
    return strategy_present
