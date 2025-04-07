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
    This function detects a synthetic strategy involving preservation of a trifluoromethyl group
    throughout the synthesis while modifying other functional groups.
    """
    has_cf3_throughout = True
    modification_count = 0

    def dfs_traverse(node):
        nonlocal has_cf3_throughout, modification_count

        if node["type"] == "mol" and node.get("in_stock", False) == False:
            # Check if molecule has CF3 group
            mol = Chem.MolFromSmiles(node["smiles"])
            cf3_pattern = Chem.MolFromSmarts("[c]-[C]([F])([F])[F]")

            if mol and not mol.HasSubstructMatch(cf3_pattern):
                has_cf3_throughout = False
                print(f"Molecule without CF3 group: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if a functional group modification occurred
            r_mol = Chem.MolFromSmiles(reactants[0])
            p_mol = Chem.MolFromSmiles(product)

            if r_mol and p_mol:
                # Patterns for functional groups
                patterns = [
                    Chem.MolFromSmarts("[C;H1](=O)"),  # aldehyde
                    Chem.MolFromSmarts("[C;H1,H2]([OH])"),  # alcohol
                    Chem.MolFromSmarts("[C](=O)[OH]"),  # carboxylic acid
                    Chem.MolFromSmarts("[c]-[#9,#17,#35,#53]"),  # aryl halide
                    Chem.MolFromSmarts("[c]-[N;H2]"),  # aniline
                    Chem.MolFromSmarts("[c]-[N+](=O)[O-]"),  # nitro
                ]

                for pattern in patterns:
                    r_has_pattern = r_mol.HasSubstructMatch(pattern)
                    p_has_pattern = p_mol.HasSubstructMatch(pattern)

                    if r_has_pattern != p_has_pattern:
                        modification_count += 1
                        print(
                            f"Detected functional group modification: {reactants[0]} -> {product}"
                        )
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if CF3 is preserved throughout and at least one modification occurred
    return has_cf3_throughout and modification_count > 0
