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
    This function detects if the synthetic route involves sulfonamide formation.
    """
    sulfonamide_formation_detected = False

    def dfs_traverse(node):
        nonlocal sulfonamide_formation_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for sulfonyl chloride in reactants
            sulfonyl_chloride_pattern = Chem.MolFromSmarts("[#6]-[S](=O)(=O)Cl")

            # Check for sulfonamide in product
            sulfonamide_pattern = Chem.MolFromSmarts("[#7]-[S](=O)(=O)-[#6]")

            # Check if any reactant has a sulfonyl chloride
            sulfonyl_chloride_present = False
            for reactant_smiles in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                if reactant_mol and reactant_mol.HasSubstructMatch(
                    sulfonyl_chloride_pattern
                ):
                    sulfonyl_chloride_present = True
                    break

            # Check if product has a sulfonamide
            product_mol = Chem.MolFromSmiles(product_smiles)
            if (
                sulfonyl_chloride_present
                and product_mol
                and product_mol.HasSubstructMatch(sulfonamide_pattern)
            ):
                sulfonamide_formation_detected = True
                print(f"Detected sulfonamide formation: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Sulfonamide formation detected: {sulfonamide_formation_detected}")
    return sulfonamide_formation_detected
