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
    This function detects if the synthesis route includes a Wittig reaction
    to transform a ketone to an alkene.
    """
    wittig_found = False

    def dfs_traverse(node):
        nonlocal wittig_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ketone pattern in reactants
            ketone_pattern = Chem.MolFromSmarts("[C](=O)[C]")

            # Check for phosphonium pattern (simplified)
            phosphonium_pattern = Chem.MolFromSmarts("[P+]")

            # Check for alkene pattern in product
            alkene_pattern = Chem.MolFromSmarts("[C]=[C]")

            reactant_has_ketone = False
            reactant_has_phosphonium = False
            product_has_alkene = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(ketone_pattern):
                    reactant_has_ketone = True
                if mol and mol.HasSubstructMatch(phosphonium_pattern):
                    reactant_has_phosphonium = True

            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(alkene_pattern):
                product_has_alkene = True

            if reactant_has_ketone and reactant_has_phosphonium and product_has_alkene:
                print("Wittig olefination detected")
                wittig_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return wittig_found
