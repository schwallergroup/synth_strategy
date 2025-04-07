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
    Detects a synthesis route that includes nitrile hydrolysis to form a carboxylic acid
    in the final step.
    """
    has_nitrile_hydrolysis = False

    # Define SMARTS patterns
    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
    carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")

    def dfs_traverse(node):
        nonlocal has_nitrile_hydrolysis

        if node["type"] == "reaction" and node.get("children", []):
            # Check if this is the first reaction (depth 0)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitrile hydrolysis
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(carboxylic_acid_pattern):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(nitrile_pattern):
                            has_nitrile_hydrolysis = True
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Nitrile to acid via hydrolysis: {has_nitrile_hydrolysis}")

    return has_nitrile_hydrolysis
