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
    Detects a synthesis that includes an ester hydrolysis step to form a carboxylic acid.
    """
    has_ester_hydrolysis = False

    def dfs_traverse(node):
        nonlocal has_ester_hydrolysis

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
            acid_pattern = Chem.MolFromSmarts("[C](=[O])[O;H]")

            has_ester = any(r is not None and r.HasSubstructMatch(ester_pattern) for r in reactants)
            has_acid = product is not None and product.HasSubstructMatch(acid_pattern)

            if has_ester and has_acid:
                has_ester_hydrolysis = True
                print("Detected ester hydrolysis")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_ester_hydrolysis
