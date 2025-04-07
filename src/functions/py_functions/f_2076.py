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
    Detects a sequence of functional group transformations from
    carboxylic acid → alcohol → aldehyde in consecutive reactions.
    """
    # Track if we've found each step in the sequence
    found_acid_to_alcohol = False
    found_alcohol_to_aldehyde = False

    def dfs_traverse(node, depth=0):
        nonlocal found_acid_to_alcohol, found_alcohol_to_aldehyde

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Patterns for functional groups
                carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
                alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
                aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")

                # Check for carboxylic acid → alcohol transformation
                reactant_has_acid = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(carboxylic_acid_pattern):
                        reactant_has_acid = True
                        break

                product_mol = Chem.MolFromSmiles(product)
                if (
                    reactant_has_acid
                    and product_mol
                    and product_mol.HasSubstructMatch(alcohol_pattern)
                ):
                    found_acid_to_alcohol = True
                    print(
                        f"Found carboxylic acid to alcohol transformation at depth {depth}"
                    )

                # Check for alcohol → aldehyde transformation
                reactant_has_alcohol = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(alcohol_pattern):
                        reactant_has_alcohol = True
                        break

                if (
                    reactant_has_alcohol
                    and product_mol
                    and product_mol.HasSubstructMatch(aldehyde_pattern)
                ):
                    found_alcohol_to_aldehyde = True
                    print(f"Found alcohol to aldehyde transformation at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    # Return True only if both transformations are found
    return found_acid_to_alcohol and found_alcohol_to_aldehyde
