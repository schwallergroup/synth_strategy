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
    This function detects if the synthesis involves coupling of aromatic and alicyclic fragments
    via olefination.
    """
    found_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_coupling

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aromatic pattern
                aromatic_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1")
                # Check for alicyclic pattern
                alicyclic_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C][C]1")
                # Check for olefin linker pattern
                olefin_linker = Chem.MolFromSmarts("[c]-[C]=[C]-[C]")

                product_mol = Chem.MolFromSmiles(product) if product else None

                if (
                    product_mol
                    and product_mol.HasSubstructMatch(aromatic_pattern)
                    and product_mol.HasSubstructMatch(alicyclic_pattern)
                    and product_mol.HasSubstructMatch(olefin_linker)
                ):
                    print(
                        "Found aromatic-alicyclic coupling via olefin at depth", depth
                    )
                    found_coupling = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return found_coupling
