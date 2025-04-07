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
    This function detects a strategy where a sulfonamide is formed from a sulfonic acid
    through a sulfonyl chloride intermediate.
    """
    # Track if we've seen each functional group transformation in the correct order
    seen_sulfonic_acid = False
    seen_sulfonyl_chloride = False
    seen_sulfonamide = False

    # SMARTS patterns for functional groups
    sulfonic_acid_pattern = Chem.MolFromSmarts("[#6][SX4](=[OX1])(=[OX1])[OX2H]")
    sulfonyl_chloride_pattern = Chem.MolFromSmarts("[#6][SX4](=[OX1])(=[OX1])[ClX1]")
    sulfonamide_pattern = Chem.MolFromSmarts("[#6][SX4](=[OX1])(=[OX1])[NX3]")

    def dfs_traverse(node):
        nonlocal seen_sulfonic_acid, seen_sulfonyl_chloride, seen_sulfonamide

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for sulfonamide in final product
                if node.get("in_stock", False) == False and mol.HasSubstructMatch(
                    sulfonamide_pattern
                ):
                    seen_sulfonamide = True
                    print("Found sulfonamide in product")

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                reactant_mol = Chem.MolFromSmiles(reactants)
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    # Check for sulfonic acid to sulfonyl chloride
                    if reactant_mol.HasSubstructMatch(
                        sulfonic_acid_pattern
                    ) and product_mol.HasSubstructMatch(sulfonyl_chloride_pattern):
                        seen_sulfonic_acid = True
                        seen_sulfonyl_chloride = True
                        print("Found sulfonic acid to sulfonyl chloride conversion")

                    # Check for sulfonyl chloride to sulfonamide
                    elif reactant_mol.HasSubstructMatch(
                        sulfonyl_chloride_pattern
                    ) and product_mol.HasSubstructMatch(sulfonamide_pattern):
                        seen_sulfonyl_chloride = True
                        seen_sulfonamide = True
                        print("Found sulfonyl chloride to sulfonamide conversion")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we've seen the complete sequence
    strategy_present = seen_sulfonic_acid and seen_sulfonyl_chloride and seen_sulfonamide
    print(f"Sulfonamide from sulfonic acid strategy detected: {strategy_present}")
    return strategy_present
