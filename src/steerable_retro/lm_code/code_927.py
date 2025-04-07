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
    This function detects if the synthesis route involves a nitro reduction to amine.
    """
    nitro_reduction_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactant contains nitro group
            nitro_pattern = Chem.MolFromSmarts("[#6][N+](=[O])[O-]")

            # Check if product contains amine group
            amine_pattern = Chem.MolFromSmarts("[#6][NH2]")

            for reactant in reactants:
                react_mol = Chem.MolFromSmiles(reactant)
                prod_mol = Chem.MolFromSmiles(product)

                if (
                    react_mol
                    and prod_mol
                    and react_mol.HasSubstructMatch(nitro_pattern)
                    and prod_mol.HasSubstructMatch(amine_pattern)
                ):
                    print(f"Nitro reduction detected at depth {depth}")
                    nitro_reduction_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return nitro_reduction_detected
