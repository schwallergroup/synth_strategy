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
    This function detects a synthetic strategy involving ester hydrolysis.
    """
    found_ester_hydrolysis = False

    def dfs_traverse(node):
        nonlocal found_ester_hydrolysis

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            products = rsmi.split(">")[-1].split(".")

            # Check for ester in reactants
            reactant_has_ester = False
            for reactant in reactants:
                if reactant:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[C]-[O]-[C](=[O])")):
                            reactant_has_ester = True
                            break
                    except:
                        continue

            # Check for alcohol in products
            product_has_alcohol = False
            for product in products:
                if product:
                    try:
                        mol = Chem.MolFromSmiles(product)
                        if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[O;H1]-[C]")):
                            product_has_alcohol = True
                            break
                    except:
                        continue

            if reactant_has_ester and product_has_alcohol:
                found_ester_hydrolysis = True
                print(f"Detected ester hydrolysis: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Ester hydrolysis strategy detected: {found_ester_hydrolysis}")
    return found_ester_hydrolysis
