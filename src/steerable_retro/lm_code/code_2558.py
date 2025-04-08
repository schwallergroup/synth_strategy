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
    This function detects an SNAr transformation strategy where a chloro group
    on a heterocycle is replaced with an amino group.
    """
    snar_found = False

    def dfs_traverse(node, depth=0):
        nonlocal snar_found

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for chloro-heterocycle to amino-heterocycle transformation
            chloro_pattern = "[c][Cl]"
            amino_pattern = "[c][NH2]"

            # Check if there's a chloro group in reactants and amino group in product
            has_chloro_reactant = any(
                Chem.MolFromSmiles(r)
                and Chem.MolFromSmiles(r).HasSubstructMatch(Chem.MolFromSmarts(chloro_pattern))
                for r in reactants
                if Chem.MolFromSmiles(r)
            )

            product_mol = Chem.MolFromSmiles(product)
            has_amino_product = product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts(amino_pattern)
            )

            if has_chloro_reactant and has_amino_product:
                print(f"Found SNAr transformation (Cl → NH2) at depth {depth}")
                snar_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return snar_found
