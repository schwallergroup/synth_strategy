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
    This function detects a strategy involving transformation of an aryl halide
    (specifically bromide) to an aryl ester.
    """
    halogen_to_ester_found = False

    def dfs_traverse(node):
        nonlocal halogen_to_ester_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl bromide to aryl ester transformation
            bromide_pattern = Chem.MolFromSmarts("[#6]:[#6]-[#35]")
            ester_pattern = Chem.MolFromSmarts("[#6]:[#6]-[#6](=[#8])-[#8]-[#6]")

            reactants_have_bromide = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(bromide_pattern)
                for r in reactants
                if Chem.MolFromSmiles(r)
            )
            product_has_ester = (
                Chem.MolFromSmiles(product).HasSubstructMatch(ester_pattern)
                if Chem.MolFromSmiles(product)
                else False
            )

            if reactants_have_bromide and product_has_ester:
                halogen_to_ester_found = True
                print("Halogen to ester transformation detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return halogen_to_ester_found
