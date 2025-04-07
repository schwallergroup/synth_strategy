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
    This function detects the use of nitrile hydrolysis to form carboxylic acids.
    """
    nitrile_hydrolysis_found = False

    def dfs_traverse(node):
        nonlocal nitrile_hydrolysis_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # Check if any reactant has a nitrile group
                reactant_has_nitrile = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[C]#[N]")):
                        reactant_has_nitrile = True
                        break

                # Check if product has a carboxylic acid group
                product_mol = Chem.MolFromSmiles(product)
                product_has_acid = product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[C](=O)[OH]")
                )

                # If reactant has nitrile and product has acid, it's likely a nitrile hydrolysis
                if reactant_has_nitrile and product_has_acid:
                    nitrile_hydrolysis_found = True
                    print(f"Nitrile hydrolysis detected in reaction: {rsmi}")
            except:
                print(f"Error processing reaction: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Nitrile hydrolysis to acid detected: {nitrile_hydrolysis_found}")
    return nitrile_hydrolysis_found
