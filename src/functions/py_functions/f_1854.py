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
    This function detects routes with ester hydrolysis.
    """
    has_ester_hydrolysis = False

    def dfs_traverse(node):
        nonlocal has_ester_hydrolysis

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for ester hydrolysis: ester in reactant, carboxylic acid in product
            for reactant_smiles in reactants_smiles:
                reactant = Chem.MolFromSmiles(reactant_smiles)
                product = Chem.MolFromSmiles(product_smiles)

                if reactant is not None and product is not None:
                    ester_pattern = Chem.MolFromSmarts("[C](=O)[O][C]")
                    carboxylic_pattern = Chem.MolFromSmarts("[C](=O)[OH]")

                    if reactant.HasSubstructMatch(
                        ester_pattern
                    ) and product.HasSubstructMatch(carboxylic_pattern):
                        has_ester_hydrolysis = True
                        print(
                            f"Detected ester hydrolysis at depth {node.get('depth', 'unknown')}"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_ester_hydrolysis
