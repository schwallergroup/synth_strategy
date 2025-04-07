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
    Detects nitro group reduction to amine as a key functional group transformation.
    """
    nitro_reduction_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro group in reactants
            nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
            has_nitro = False
            for r in reactants:
                r_mol = Chem.MolFromSmiles(r)
                if r_mol and r_mol.HasSubstructMatch(nitro_pattern):
                    has_nitro = True
                    break

            # Check for amine in product where nitro was
            if has_nitro:
                product_mol = Chem.MolFromSmiles(product)
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                    # This is a simplification - in a real implementation,
                    # we would need to verify the amine is at the same position as the nitro was
                    nitro_reduction_detected = True
                    print(f"Detected nitro reduction to amine at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return nitro_reduction_detected
