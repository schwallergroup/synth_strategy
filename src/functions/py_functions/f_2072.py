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
    This function detects if the synthesis involves incorporation of a sulfonamide group.
    """
    sulfonamide_incorporated = False

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_incorporated

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                # Check for sulfonamide pattern in reactants
                sulfonamide_pattern = Chem.MolFromSmarts("[#7]-[S](=[O])(=[O])-[#6]")

                # Check if sulfonamide is present in product but not in all reactants
                if product_mol and product_mol.HasSubstructMatch(sulfonamide_pattern):
                    for r_mol in reactant_mols:
                        if r_mol and r_mol.HasSubstructMatch(sulfonamide_pattern):
                            print(f"Found sulfonamide incorporation at depth {depth}")
                            sulfonamide_incorporated = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return sulfonamide_incorporated
