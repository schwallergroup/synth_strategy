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
    This function detects if the synthetic route contains a sulfonamide formation step.
    """
    sulfonamide_formation_found = False

    def dfs_traverse(node):
        nonlocal sulfonamide_formation_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if rsmi:
                product_smiles = rsmi.split(">")[-1]

                # Check if product has sulfonamide group
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    sulfonamide_pattern = Chem.MolFromSmarts("[NH]-[S](=O)(=O)")

                    # Check if this reaction forms the sulfonamide
                    if product_mol.HasSubstructMatch(sulfonamide_pattern):
                        reactants_smiles = rsmi.split(">")[0].split(".")

                        # Check if reactants don't have the complete sulfonamide pattern
                        sulfonamide_in_reactants = False
                        for r in reactants_smiles:
                            r_mol = Chem.MolFromSmiles(r)
                            if r_mol and r_mol.HasSubstructMatch(sulfonamide_pattern):
                                sulfonamide_in_reactants = True

                        if not sulfonamide_in_reactants:
                            print("Sulfonamide formation detected")
                            sulfonamide_formation_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return sulfonamide_formation_found
