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
    This function detects if the synthesis uses Boc protection strategy
    for amine groups.
    """
    boc_protection_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal boc_protection_count

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a Boc protection
                product_mol = Chem.MolFromSmiles(product)
                boc_pattern = Chem.MolFromSmarts(
                    "[#7]-[#6](=[#8])-[#8]-[#6]([#6])([#6])[#6]"
                )

                if product_mol and product_mol.HasSubstructMatch(boc_pattern):
                    # Check if this is a new Boc group (not present in all reactants)
                    has_new_boc = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and not reactant_mol.HasSubstructMatch(
                            boc_pattern
                        ):
                            has_new_boc = True
                            break

                    if has_new_boc:
                        print(f"Found Boc protection at depth {depth}")
                        boc_protection_count += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return boc_protection_count >= 2  # At least 2 Boc protections
