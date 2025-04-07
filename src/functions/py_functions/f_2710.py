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
    Detects if the synthesis route involves a Suzuki coupling reaction.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for boronic acid in reactants
            has_boronic_acid = any(
                "B(O)(O)" in r or "OB(O)" in r for r in reactants_smiles
            )

            # Check for halogen (Cl, Br, I) in reactants
            has_halogen = any(
                any(x in r for x in ["Cl", "Br", "I"]) for r in reactants_smiles
            )

            # Check for aromatic rings in reactants and product
            has_aromatic_reactants = any("c" in r or "n" in r for r in reactants_smiles)
            has_aromatic_product = "c" in product_smiles or "n" in product_smiles

            if (
                has_boronic_acid
                and has_halogen
                and has_aromatic_reactants
                and has_aromatic_product
            ):
                print("Detected Suzuki coupling reaction")
                result = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return result
