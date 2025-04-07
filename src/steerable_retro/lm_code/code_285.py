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
    This function detects the use of sulfonylmethyl (CH2-SO2) as a linking group between aromatic fragments.
    """
    found_sulfonylmethyl_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_sulfonylmethyl_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for sulfonylmethyl linker in product
            sulfonylmethyl_pattern = "[c][CH2][S](=O)(=O)[c]"

            product = Chem.MolFromSmiles(product_smiles)
            if product and product.HasSubstructMatch(Chem.MolFromSmarts(sulfonylmethyl_pattern)):
                # Check if this pattern is new (not in all reactants)
                all_have_pattern = True
                for r_smiles in reactants_smiles:
                    r = Chem.MolFromSmiles(r_smiles)
                    if r and not r.HasSubstructMatch(Chem.MolFromSmarts(sulfonylmethyl_pattern)):
                        all_have_pattern = False
                        break

                if not all_have_pattern:
                    found_sulfonylmethyl_formation = True
                    print(f"Sulfonylmethyl linker formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_sulfonylmethyl_formation
