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
    Detects if the synthesis route involves parallel transformations of nitrile groups
    to different functional groups (carboxylic acid and oxime).
    """
    nitrile_to_acid = False
    nitrile_to_oxime = False

    def dfs_traverse(node):
        nonlocal nitrile_to_acid, nitrile_to_oxime

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitrile to carboxylic acid transformation
                nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(acid_pattern):
                    for r_mol in reactant_mols:
                        if r_mol and r_mol.HasSubstructMatch(nitrile_pattern):
                            nitrile_to_acid = True
                            print("Detected nitrile to carboxylic acid transformation")

                # Check for nitrile to oxime transformation
                oxime_pattern = Chem.MolFromSmarts("[C](=[N][OH])")

                if product_mol and product_mol.HasSubstructMatch(oxime_pattern):
                    for r_mol in reactant_mols:
                        if r_mol and r_mol.HasSubstructMatch(nitrile_pattern):
                            nitrile_to_oxime = True
                            print("Detected nitrile to oxime transformation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both transformations are found
    return nitrile_to_acid and nitrile_to_oxime
