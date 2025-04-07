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
    Detects if the synthesis route involves sequential redox transformations,
    specifically looking for carboxylic acid → alcohol → aldehyde sequence.
    """
    # Track if we've seen each transformation
    acid_to_alcohol = False
    alcohol_to_aldehyde = False

    def dfs_traverse(node, depth=0):
        nonlocal acid_to_alcohol, alcohol_to_aldehyde

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Define patterns
            carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
            alcohol_pattern = Chem.MolFromSmarts("[CX4H2][OX2H]")
            aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=[OX1])[#6]")

            # Check for acid to alcohol transformation
            try:
                reactant_mol = Chem.MolFromSmiles(reactants[0])
                product_mol = Chem.MolFromSmiles(product)

                if (
                    reactant_mol
                    and product_mol
                    and reactant_mol.HasSubstructMatch(carboxylic_acid_pattern)
                    and product_mol.HasSubstructMatch(alcohol_pattern)
                    and not reactant_mol.HasSubstructMatch(alcohol_pattern)
                ):
                    print(
                        f"Detected carboxylic acid to alcohol reduction at depth {depth}"
                    )
                    acid_to_alcohol = True

                # Check for alcohol to aldehyde transformation
                if (
                    reactant_mol
                    and product_mol
                    and reactant_mol.HasSubstructMatch(alcohol_pattern)
                    and product_mol.HasSubstructMatch(aldehyde_pattern)
                    and not reactant_mol.HasSubstructMatch(aldehyde_pattern)
                ):
                    print(f"Detected alcohol to aldehyde oxidation at depth {depth}")
                    alcohol_to_aldehyde = True
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    # Return True only if both transformations are detected
    return acid_to_alcohol and alcohol_to_aldehyde
