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
    Detects if the synthesis route involves benzimidazole ring formation.
    """
    benzimidazole_formed = False

    def dfs_traverse(node):
        nonlocal benzimidazole_formed

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for benzimidazole formation
                benzimidazole_pattern = Chem.MolFromSmarts("[nH]1cnc2ccccc12")

                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(benzimidazole_pattern):
                    # Check if benzimidazole is not in reactants
                    benzimidazole_in_reactants = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            benzimidazole_pattern
                        ):
                            benzimidazole_in_reactants = True
                            break

                    if not benzimidazole_in_reactants:
                        benzimidazole_formed = True
                        print("Detected benzimidazole formation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return benzimidazole_formed
