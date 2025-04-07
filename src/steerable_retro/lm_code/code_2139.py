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
    Detects a synthesis involving SNAr connection of heteroaryl fragments.
    """
    has_snar_connection = False

    def dfs_traverse(node):
        nonlocal has_snar_connection

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            if len(reactants) >= 2:
                # Check for phenol and chloroheteroaryl patterns in reactants
                phenol_present = False
                chloroheteroaryl_present = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        phenol_pattern = Chem.MolFromSmarts("[OH][c]")
                        if mol.HasSubstructMatch(phenol_pattern):
                            phenol_present = True

                        # Pattern for chloroheteroaryl (focusing on benzothiazole)
                        chloroheteroaryl_pattern = Chem.MolFromSmarts(
                            "[Cl][c]1[c]2[s][c][cH][c]2[n][cH][cH]1"
                        )
                        if mol.HasSubstructMatch(chloroheteroaryl_pattern):
                            chloroheteroaryl_present = True

                # Check for connected structure in product
                if phenol_present and chloroheteroaryl_present:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        connected_pattern = Chem.MolFromSmarts(
                            "[c][O][c]1[c]2[s][c][cH][c]2[n][cH][cH]1"
                        )
                        if product_mol.HasSubstructMatch(connected_pattern):
                            has_snar_connection = True
                            print(f"Found SNAr heteroaryl connection in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"SNAr heteroaryl connection strategy detected: {has_snar_connection}")
    return has_snar_connection
