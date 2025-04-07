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
    This function detects if the final step in the synthesis involves formation of a sulfonamide bond
    between a sulfonyl chloride and an amine.
    """
    final_reaction_found = False
    sulfonamide_formation = False

    def dfs_traverse(node):
        nonlocal final_reaction_found, sulfonamide_formation

        if node["type"] == "reaction" and not final_reaction_found:
            # This is the first reaction node we encounter (depth 0)
            final_reaction_found = True

            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for sulfonyl chloride pattern in reactants
            sulfonyl_chloride_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])[Cl]")
            amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")

            # Check for sulfonamide pattern in product
            sulfonamide_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])[#7]")

            # Check if reactants contain sulfonyl chloride and amine
            has_sulfonyl_chloride = False
            has_amine = False

            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(sulfonyl_chloride_pattern):
                        has_sulfonyl_chloride = True
                    if mol.HasSubstructMatch(amine_pattern):
                        has_amine = True

            # Check if product contains sulfonamide
            product_mol = Chem.MolFromSmiles(product_smiles)
            has_sulfonamide = False
            if product_mol and product_mol.HasSubstructMatch(sulfonamide_pattern):
                has_sulfonamide = True

            # If all conditions are met, this is a sulfonamide formation
            if has_sulfonyl_chloride and has_amine and has_sulfonamide:
                sulfonamide_formation = True
                print("Detected late-stage sulfonamide formation")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return sulfonamide_formation
