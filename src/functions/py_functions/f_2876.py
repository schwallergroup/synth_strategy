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
    Detects if stereocenters are established early in the synthesis (at deeper levels of the tree).
    """
    stereocenters_established = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal stereocenters_established, max_depth

        max_depth = max(max_depth, depth)

        if (
            node["type"] == "reaction" and depth >= 3
        ):  # Early in synthesis (deeper in tree)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check for stereocenters in product
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    Chem.AssignStereochemistry(product_mol)
                    chiral_centers = [
                        atom.GetIdx()
                        for atom in product_mol.GetAtoms()
                        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
                    ]

                    if len(chiral_centers) > 0:
                        print(f"Stereocenters established at depth {depth}")
                        stereocenters_established = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Only consider it early if we have stereocenters and sufficient depth
    is_early = stereocenters_established and max_depth >= 4
    print(f"Early stereocenter establishment: {is_early}")
    return is_early
