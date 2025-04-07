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
    Detects if the synthesis includes protection of an alcohol as a vinyl ether.
    """
    found_vinyl_ether_protection = False

    def dfs_traverse(node):
        nonlocal found_vinyl_ether_protection

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for vinyl ether formation
                alcohol_pattern = Chem.MolFromSmarts("[CH2][OH]")
                vinyl_ether_pattern = Chem.MolFromSmarts("[CH2][O][CH]=[CH2]")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(alcohol_pattern):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(
                            vinyl_ether_pattern
                        ):
                            print("Found vinyl ether protection of alcohol")
                            found_vinyl_ether_protection = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_vinyl_ether_protection
