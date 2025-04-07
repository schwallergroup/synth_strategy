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
    This function detects if the synthesis route involves reduction of a nitro group to an amine.
    """
    nitro_reduction_detected = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitro group in reactants
            nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            # Look for nitro group in reactants
            nitro_in_reactants = False
            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(nitro_pattern):
                    nitro_in_reactants = True

                    # Check if the product has an amine where the nitro group was
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                        # This is a simplification - ideally we would check if the amine is at the same position
                        nitro_reduction_detected = True
                        print("Nitro reduction to amine detected")
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return nitro_reduction_detected
