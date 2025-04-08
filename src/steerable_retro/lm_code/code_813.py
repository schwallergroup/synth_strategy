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
    This function detects if the synthesis involves a Wittig-type reaction
    (using a phosphonium reagent to form a C=C bond).
    """
    wittig_reaction_found = False

    def dfs_traverse(node):
        nonlocal wittig_reaction_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for phosphonium in reactants
            phosphonium_pattern = Chem.MolFromSmarts("[P+]")
            phosphonium_found = False

            # Check for new C=C in product that wasn't in reactants
            alkene_pattern = Chem.MolFromSmarts("[C]=[C]")

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(phosphonium_pattern):
                        phosphonium_found = True
                        print("Phosphonium reagent found in reactants")
                except:
                    continue

            if phosphonium_found:
                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(alkene_pattern):
                        wittig_reaction_found = True
                        print("Wittig-type reaction detected")
                except:
                    pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return wittig_reaction_found
