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
    Detects if the synthesis route involves the formation of a heterocyclic ring
    (specifically looking for pyrazine formation).
    """
    heterocycle_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for pyrazine formation
            pyrazine_pattern = Chem.MolFromSmarts("[n]1[c][n][c][c][c]1")

            # Check if product contains pyrazine
            try:
                prod_mol = Chem.MolFromSmiles(product)
                has_pyrazine = prod_mol and prod_mol.HasSubstructMatch(pyrazine_pattern)
            except:
                has_pyrazine = False

            # Check if reactants don't contain pyrazine
            reactants_have_pyrazine = False
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(pyrazine_pattern):
                        reactants_have_pyrazine = True
                        break
                except:
                    continue

            # If product has pyrazine but reactants don't, it's a pyrazine formation
            if has_pyrazine and not reactants_have_pyrazine:
                print(f"Detected heterocycle (pyrazine) formation at depth {depth}")
                heterocycle_formation = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return heterocycle_formation
