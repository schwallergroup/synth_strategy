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
    This function detects if the synthesis route involves multiple ether formation reactions,
    particularly benzyl ether formations.
    """
    ether_formation_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal ether_formation_count

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for hydroxyl group in reactants
                hydroxyl_pattern = Chem.MolFromSmarts("[#8;H1]-[#6]")
                benzyl_pattern = Chem.MolFromSmarts("[#6]-[#6]~[c]")
                ether_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]~[c]")

                has_hydroxyl = False
                has_benzyl = False

                for r in reactants:
                    if r:
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            if mol.HasSubstructMatch(hydroxyl_pattern):
                                has_hydroxyl = True
                            if mol.HasSubstructMatch(benzyl_pattern):
                                has_benzyl = True

                # Check for benzyl ether in product
                product_mol = Chem.MolFromSmiles(product)
                if (
                    product_mol
                    and has_hydroxyl
                    and has_benzyl
                    and product_mol.HasSubstructMatch(ether_pattern)
                ):
                    print(f"Ether formation detected at depth {depth}")
                    ether_formation_count += 1

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return ether_formation_count >= 2
