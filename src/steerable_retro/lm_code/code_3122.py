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
    This function detects if the synthesis involves formation of a carboxylic acid
    from an aromatic C-H in the early stages of the synthesis.
    """
    found_carboxylic_acid_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_carboxylic_acid_formation

        if node["type"] == "reaction" and depth >= 3:  # Early stage reaction
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for aromatic to carboxylic acid transformation
            reactant_is_aromatic = False
            product_has_carboxylic_acid = False

            for reactant in reactants_smiles:
                if reactant:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.GetNumAtoms() > 6:  # Ensure it's not just a small fragment
                            # Check if molecule has aromatic atoms
                            for atom in mol.GetAtoms():
                                if atom.GetIsAromatic():
                                    reactant_is_aromatic = True
                                    break
                    except:
                        continue

            if product_smiles:
                try:
                    mol = Chem.MolFromSmiles(product_smiles)
                    if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[OH]")):
                        product_has_carboxylic_acid = True
                except:
                    pass

            if reactant_is_aromatic and product_has_carboxylic_acid:
                found_carboxylic_acid_formation = True
                print(f"Found carboxylic acid formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_carboxylic_acid_formation
