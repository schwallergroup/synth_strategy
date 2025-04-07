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
    Detects if the synthesis involves formation of an aryl ether (C-O-C where one C is aromatic).
    Looks for reactions where a hydroxyl group and aryl halide form an ether.
    """
    found_ether_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_ether_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Patterns for hydroxyl and aryl halide
                hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
                aryl_halide_pattern = Chem.MolFromSmarts("c[Br,I,Cl]")

                # Check for hydroxyl and aryl halide in reactants
                has_hydroxyl = False
                has_aryl_halide = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if not mol:
                        continue
                    if mol.HasSubstructMatch(hydroxyl_pattern):
                        has_hydroxyl = True
                    if mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True

                # Check for aryl ether in product
                product_mol = Chem.MolFromSmiles(product)
                aryl_ether_pattern = Chem.MolFromSmarts("c-[#8]-[#6]")
                has_aryl_ether = product_mol and product_mol.HasSubstructMatch(aryl_ether_pattern)

                if has_hydroxyl and has_aryl_halide and has_aryl_ether:
                    found_ether_formation = True
                    print(f"Detected aryl ether formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_ether_formation
