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
    This function detects a synthetic strategy involving both oxidation and reduction reactions
    in the same route, specifically looking for carbonyl reductions and alcohol oxidations.
    """
    oxidation_count = 0
    reduction_count = 0

    def dfs_traverse(node):
        nonlocal oxidation_count, reduction_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for carbonyl reduction (aldehyde/ketone to alcohol)
            aldehyde_pattern = Chem.MolFromSmarts("[C;H1](=O)")
            ketone_pattern = Chem.MolFromSmarts("[C](=O)[C,c]")
            alcohol_pattern = Chem.MolFromSmarts("[C;H1,H2]([OH])")

            for reactant in reactants:
                r_mol = Chem.MolFromSmiles(reactant)
                if r_mol:
                    if r_mol.HasSubstructMatch(aldehyde_pattern) or r_mol.HasSubstructMatch(
                        ketone_pattern
                    ):
                        p_mol = Chem.MolFromSmiles(product)
                        if p_mol and p_mol.HasSubstructMatch(alcohol_pattern):
                            reduction_count += 1
                            print(f"Detected carbonyl reduction: {reactant} -> {product}")

            # Check for alcohol oxidation (alcohol to carboxylic acid or aldehyde)
            alcohol_pattern = Chem.MolFromSmarts("[C;H1,H2]([OH])")
            carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")
            aldehyde_pattern = Chem.MolFromSmarts("[C;H1](=O)")

            for reactant in reactants:
                r_mol = Chem.MolFromSmiles(reactant)
                if r_mol and r_mol.HasSubstructMatch(alcohol_pattern):
                    p_mol = Chem.MolFromSmiles(product)
                    if p_mol and (
                        p_mol.HasSubstructMatch(carboxylic_acid_pattern)
                        or p_mol.HasSubstructMatch(aldehyde_pattern)
                    ):
                        oxidation_count += 1
                        print(f"Detected alcohol oxidation: {reactant} -> {product}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both oxidation and reduction are present
    return oxidation_count > 0 and reduction_count > 0
