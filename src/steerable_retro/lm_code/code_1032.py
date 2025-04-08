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
    Detects if the synthesis includes an ester hydrolysis step.
    """
    found_ester_hydrolysis = False

    def dfs_traverse(node):
        nonlocal found_ester_hydrolysis

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for ester hydrolysis
            ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]")
            acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")

            if any(
                r.HasSubstructMatch(ester_pattern) for r in reactants
            ) and product.HasSubstructMatch(acid_pattern):
                found_ester_hydrolysis = True
                print("Found ester hydrolysis")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_ester_hydrolysis
