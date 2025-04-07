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
    This function detects if the synthesis follows a linear fragment combination strategy
    with at least 2 fragment combinations.
    """
    fragment_combinations = 0

    def dfs_traverse(node):
        nonlocal fragment_combinations

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Count number of fragments in reactants and product
            if len(reactants_smiles) > 1:  # More than one reactant
                try:
                    p_mol = Chem.MolFromSmiles(product_smiles)
                    if p_mol:
                        # If multiple reactants combine to form a single product, it's a fragment combination
                        fragment_combinations += 1
                        print(
                            f"Detected fragment combination (count: {fragment_combinations})"
                        )
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return fragment_combinations >= 2  # Return True if at least 2 fragment combinations
