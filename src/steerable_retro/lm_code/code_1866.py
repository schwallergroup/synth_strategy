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
    This function detects N-oxide formation and reduction in pyridine synthesis.
    Looks for pyridine N-oxide intermediates.
    """
    n_oxide_formation = False
    n_oxide_reduction = False

    def dfs_traverse(node):
        nonlocal n_oxide_formation, n_oxide_reduction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # N-oxide pattern
            n_oxide_pattern = Chem.MolFromSmarts("[n+]([O-])")
            pyridine_pattern = Chem.MolFromSmarts("n1ccccc1")

            # Check for N-oxide formation
            product_mol = Chem.MolFromSmiles(product)
            reactants_with_n_oxide = 0

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(n_oxide_pattern):
                    reactants_with_n_oxide += 1

            if (
                product_mol
                and product_mol.HasSubstructMatch(n_oxide_pattern)
                and reactants_with_n_oxide == 0
            ):
                n_oxide_formation = True
                print(f"N-oxide formation detected in product: {product}")

            # Check for N-oxide reduction
            reactants_with_n_oxide = 0
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(n_oxide_pattern):
                    reactants_with_n_oxide += 1

            if (
                reactants_with_n_oxide > 0
                and product_mol
                and product_mol.HasSubstructMatch(pyridine_pattern)
                and not product_mol.HasSubstructMatch(n_oxide_pattern)
            ):
                n_oxide_reduction = True
                print(f"N-oxide reduction detected: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return n_oxide_formation and n_oxide_reduction
