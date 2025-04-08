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
    This function detects if the synthesis involves formation of a nitrile from an amide,
    which itself was formed from an acid chloride.
    """
    found_amide_to_nitrile = False
    found_acid_chloride_to_amide = False

    def dfs_traverse(node):
        nonlocal found_amide_to_nitrile, found_acid_chloride_to_amide

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for amide to nitrile transformation
            reactant_has_amide = False
            product_has_nitrile = False

            for reactant in reactants_smiles:
                if reactant:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[NX3]")):
                            reactant_has_amide = True
                            break
                    except:
                        continue

            if product_smiles:
                try:
                    mol = Chem.MolFromSmiles(product_smiles)
                    if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[CX2]#[NX1]")):
                        product_has_nitrile = True
                except:
                    pass

            if reactant_has_amide and product_has_nitrile:
                found_amide_to_nitrile = True
                print("Found amide to nitrile transformation")

            # Check for acid chloride to amide transformation
            reactant_has_acid_chloride = False
            product_has_amide = False

            for reactant in reactants_smiles:
                if reactant:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[Cl]")):
                            reactant_has_acid_chloride = True
                            break
                    except:
                        continue

            if product_smiles:
                try:
                    mol = Chem.MolFromSmiles(product_smiles)
                    if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[NX3]")):
                        product_has_amide = True
                except:
                    pass

            if reactant_has_acid_chloride and product_has_amide:
                found_acid_chloride_to_amide = True
                print("Found acid chloride to amide transformation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_amide_to_nitrile and found_acid_chloride_to_amide
