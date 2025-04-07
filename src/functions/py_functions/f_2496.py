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
    Detects if the synthesis route involves early-stage functional group manipulations
    (at depths 4-5) before heterocycle formation and coupling.
    """
    early_fg_manipulation = False

    def dfs_traverse(node, depth=0):
        nonlocal early_fg_manipulation

        if node["type"] == "reaction" and depth >= 4:
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Define patterns for functional group transformations
            carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
            alcohol_pattern = Chem.MolFromSmarts("[CX4H2][OX2H]")
            aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=[OX1])[#6]")

            try:
                reactant_mol = Chem.MolFromSmiles(reactants[0])
                product_mol = Chem.MolFromSmiles(product)

                # Check for acid to alcohol or alcohol to aldehyde transformations
                if reactant_mol and product_mol:
                    if reactant_mol.HasSubstructMatch(
                        carboxylic_acid_pattern
                    ) and product_mol.HasSubstructMatch(alcohol_pattern):
                        print(
                            f"Detected early carboxylic acid to alcohol reduction at depth {depth}"
                        )
                        early_fg_manipulation = True

                    if reactant_mol.HasSubstructMatch(
                        alcohol_pattern
                    ) and product_mol.HasSubstructMatch(aldehyde_pattern):
                        print(
                            f"Detected early alcohol to aldehyde oxidation at depth {depth}"
                        )
                        early_fg_manipulation = True
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return early_fg_manipulation
