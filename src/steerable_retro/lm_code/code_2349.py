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
    This function detects a synthetic strategy involving azide reduction to amine.
    It looks for a reaction where an azide group is converted to an amine.
    """
    azide_to_amine_found = False

    def dfs_traverse(node):
        nonlocal azide_to_amine_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for azide in reactants
                azide_pattern = Chem.MolFromSmarts("[N-]=[N+]=N")
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(azide_pattern):
                        # Check for amine in product
                        amine_pattern = Chem.MolFromSmarts("[NH2]C")
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                            print("Found azide to amine conversion")
                            azide_to_amine_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return azide_to_amine_found
