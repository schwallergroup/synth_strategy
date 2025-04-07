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
    This function detects a synthetic strategy involving sequential borylation reactions
    (aryl halide to boronic acid/ester transformations).
    """
    borylation_count = 0

    def dfs_traverse(node):
        nonlocal borylation_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for borylation: aryl halide to boronic acid/ester
            reactant_has_aryl_halide = False
            product_has_boronic = False

            for reactant in reactants:
                if reactant:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[Br,I,Cl]")):
                            reactant_has_aryl_halide = True
                            break
                    except:
                        continue

            if product:
                try:
                    mol = Chem.MolFromSmiles(product)
                    if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[B]([O])([O])")):
                        product_has_boronic = True
                except:
                    pass

            if reactant_has_aryl_halide and product_has_boronic:
                borylation_count += 1
                print(f"Detected borylation reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is present if we have at least 2 borylation reactions
    result = borylation_count >= 2
    print(f"Borylation sequence strategy detected: {result} (count: {borylation_count})")
    return result
