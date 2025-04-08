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
    This function detects if the synthesis involves benzylic functionalization,
    particularly benzylic bromination.
    """
    benzylic_functionalization_detected = False

    def dfs_traverse(node):
        nonlocal benzylic_functionalization_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for benzylic methyl to bromomethyl transformation
                methyl_pattern = Chem.MolFromSmarts("[c][CH3]")
                bromomethyl_pattern = Chem.MolFromSmarts("[c][CH2][Br]")

                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(bromomethyl_pattern):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(methyl_pattern):
                            print("Benzylic bromination detected")
                            benzylic_functionalization_detected = True
                            break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    if benzylic_functionalization_detected:
        print("Benzylic functionalization strategy detected")
    else:
        print("No benzylic functionalization strategy detected")

    return benzylic_functionalization_detected
