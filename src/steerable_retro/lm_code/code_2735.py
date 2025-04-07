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
    This function detects a synthetic strategy involving activation of carboxylic acid
    as acyl chloride for subsequent coupling.
    """
    # Track if we found the acyl chloride formation
    found_acyl_chloride_formation = False

    def dfs_traverse(node):
        nonlocal found_acyl_chloride_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for carboxylic acid to acyl chloride conversion
            carboxylic_acid_pattern = Chem.MolFromSmarts("[C$(C=O)][OH]")
            acyl_chloride_pattern = Chem.MolFromSmarts("[C$(C=O)][Cl]")

            # Check reactants for carboxylic acid
            has_carboxylic_acid = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(carboxylic_acid_pattern):
                    has_carboxylic_acid = True
                    break

            # Check product for acyl chloride
            product_mol = Chem.MolFromSmiles(product)
            has_acyl_chloride = product_mol and product_mol.HasSubstructMatch(acyl_chloride_pattern)

            if has_carboxylic_acid and has_acyl_chloride:
                found_acyl_chloride_formation = True
                print("Found acyl chloride activation of carboxylic acid")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_acyl_chloride_formation
