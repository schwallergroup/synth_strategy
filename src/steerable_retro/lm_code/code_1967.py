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
    This function detects a synthetic strategy involving N-Boc protection
    that is maintained through multiple steps.
    """
    boc_protection_steps = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains Boc group
                product_mol = Chem.MolFromSmiles(product)
                boc_pattern = Chem.MolFromSmarts("[#7]-C(=O)OC(C)(C)C")

                if product_mol and product_mol.HasSubstructMatch(boc_pattern):
                    # Check if any reactant has the Boc group
                    has_boc_in_reactants = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(boc_pattern):
                            has_boc_in_reactants = True
                            break

                    if not has_boc_in_reactants:
                        # This is a Boc protection step
                        boc_protection_steps.append(depth)
                        print(f"N-Boc protection detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if Boc protection occurs early in the synthesis (higher depth)
    # and is maintained for multiple steps
    if boc_protection_steps and max(boc_protection_steps) >= 4:
        print(f"N-Boc protection strategy detected at depths: {boc_protection_steps}")
        return True

    return False
