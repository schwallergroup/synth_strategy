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
    Detects a sequence of functional group activations:
    ester hydrolysis followed by acid chloride formation for amide coupling
    """
    # Initialize flags
    has_ester_hydrolysis = False
    has_acid_chloride_formation = False
    has_amide_formation = False

    def dfs_traverse(node):
        nonlocal has_ester_hydrolysis, has_acid_chloride_formation, has_amide_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ester hydrolysis
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[C](=[O])[O][C]")
                ):
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[C](=[O])[O]")
                    ):
                        print("Found ester hydrolysis")
                        has_ester_hydrolysis = True

            # Check for acid chloride formation
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[C](=[O])[O]")
                ):
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[C](=[O])[Cl]")
                    ):
                        print("Found acid chloride formation")
                        has_acid_chloride_formation = True

            # Check for amide formation
            if len(reactants) >= 2:
                has_acid_chloride = False
                has_amine = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[Cl]")):
                            has_acid_chloride = True
                        if mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]")):
                            has_amine = True

                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[C](=[O])[NH]")
                ):
                    if has_acid_chloride and has_amine:
                        print("Found amide formation")
                        has_amide_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if the sequence is found
    return has_ester_hydrolysis and has_acid_chloride_formation and has_amide_formation
