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
    This function detects if the synthesis includes borylation of an aryl halide.
    """
    has_borylation = False

    def dfs_traverse(node):
        nonlocal has_borylation

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl halide in reactant
            aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Br,I]")
            boronic_ester_pattern = Chem.MolFromSmarts("[c]-[B]([O])[O]")

            has_aryl_halide = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                    has_aryl_halide = True
                    break

            # Check if product has boronic ester
            if has_aryl_halide:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(boronic_ester_pattern):
                    print("Detected borylation of aryl halide")
                    has_borylation = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return has_borylation
