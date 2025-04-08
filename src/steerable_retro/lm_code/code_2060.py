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
    This function detects if a synthetic route includes transformation of aryl halides to aryl ethers.
    """
    has_transformation = False

    def dfs_traverse(node):
        nonlocal has_transformation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for aryl halide to aryl ether transformation
            reactant_mol = Chem.MolFromSmiles(reactants)
            product_mol = Chem.MolFromSmiles(product)

            if reactant_mol and product_mol:
                # Pattern for aryl halide in reactant
                aryl_halide_patt = Chem.MolFromSmarts("[c]-[#53,#35,#17,#9]")
                # Pattern for aryl ether in product
                aryl_ether_patt = Chem.MolFromSmarts("[c]-[#8]-[#6]")

                if reactant_mol.HasSubstructMatch(
                    aryl_halide_patt
                ) and product_mol.HasSubstructMatch(aryl_ether_patt):
                    print(f"Found aryl halide to ether transformation: {rsmi}")
                    has_transformation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    print(f"Aryl halide to ether transformation strategy detected: {has_transformation}")
    return has_transformation
