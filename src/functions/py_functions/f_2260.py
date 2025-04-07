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
    Detects if the synthesis involves the formation of an amide bond between
    a carboxylic acid and an amine.
    """
    # Track if we found amide formation
    found_amide_formation = False

    def is_amide_formation(reaction_smiles):
        """Check if a reaction is an amide formation"""
        # Split reaction into reactants and products
        parts = reaction_smiles.split(">")
        if len(parts) != 3:
            return False

        reactants = parts[0].split(".")
        product = parts[2]

        # Check for carboxylic acid and amine in reactants
        has_carboxylic_acid = False
        has_amine = False

        for reactant in reactants:
            try:
                mol = Chem.MolFromSmiles(reactant)
                if not mol:
                    continue

                if mol.HasSubstructMatch(Chem.MolFromSmarts("[O]-[C](=[O])")):
                    has_carboxylic_acid = True
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[N]")):
                    has_amine = True
            except:
                continue

        # Check for amide in product
        try:
            prod_mol = Chem.MolFromSmiles(product)
            has_amide = prod_mol and prod_mol.HasSubstructMatch(
                Chem.MolFromSmarts("[N]-[C](=[O])")
            )

            # If reactants have carboxylic acid and amine, and product has amide, likely an amide formation
            if has_carboxylic_acid and has_amine and has_amide:
                return True
        except:
            pass

        return False

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_formation

        if node["type"] == "reaction":
            # Get reaction SMILES
            if "rsmi" in node.get("metadata", {}):
                reaction_smiles = node["metadata"]["rsmi"]

                # Check if this is an amide formation
                if is_amide_formation(reaction_smiles):
                    found_amide_formation = True
                    print(f"Found amide formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found amide formation
    if found_amide_formation:
        print("Detected amide bond formation strategy")
        return True

    return False
