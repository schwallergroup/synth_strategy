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
    This function detects a linear synthesis strategy that combines multiple fragments
    through sequential heteroatom linkages (O, N, S).
    """
    heteroatom_linkages = []

    def dfs_traverse(node):
        nonlocal heteroatom_linkages

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for various heteroatom linkage formations
            diaryl_ether_pattern = Chem.MolFromSmarts("c-[O]-c")
            amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
            sulfonamide_pattern = Chem.MolFromSmarts("[S](=[O])(=[O])[N]")

            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol:
                if product_mol.HasSubstructMatch(diaryl_ether_pattern):
                    heteroatom_linkages.append("O")
                    print("Found diaryl ether linkage (O)")

                if product_mol.HasSubstructMatch(amide_pattern):
                    heteroatom_linkages.append("N")
                    print("Found amide linkage (N)")

                if product_mol.HasSubstructMatch(sulfonamide_pattern):
                    heteroatom_linkages.append("S")
                    print("Found sulfonamide linkage (S)")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Check if we have at least 3 different heteroatom linkages in a linear sequence
    unique_linkages = set(heteroatom_linkages)
    return len(unique_linkages) >= 3
