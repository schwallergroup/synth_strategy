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
    Detects a synthesis route that includes a nitro reduction to amine step.
    """
    has_nitro_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal has_nitro_reduction

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check for nitro reduction
            nitro_pattern = Chem.MolFromSmarts("[#7+](=[#8])[#8-]")
            amine_pattern = Chem.MolFromSmarts("[#7H2]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
            product_mol = Chem.MolFromSmiles(product_part)

            if (
                product_mol
                and product_mol.HasSubstructMatch(amine_pattern)
                and any(r and r.HasSubstructMatch(nitro_pattern) for r in reactant_mols if r)
            ):
                print(f"Found nitro reduction at depth {depth}")
                has_nitro_reduction = True

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    if has_nitro_reduction:
        print("Nitro reduction to amine sequence detected")
    else:
        print("Nitro reduction to amine sequence NOT detected")

    return has_nitro_reduction
