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
    Detects if the synthesis includes activation of carboxylic acid to acid chloride.
    """
    found_acid_activation = False

    def dfs_traverse(node):
        nonlocal found_acid_activation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acid chloride formation pattern
                acid_chloride_pattern = Chem.MolFromSmarts("[#6:1](=[O:2])[Cl:3]")
                prod_mol = Chem.MolFromSmiles(product)

                if prod_mol and prod_mol.HasSubstructMatch(acid_chloride_pattern):
                    for reactant in reactants:
                        carboxylic_acid_pattern = Chem.MolFromSmarts("[#6](=[O])[OH]")
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol and r_mol.HasSubstructMatch(carboxylic_acid_pattern):
                            found_acid_activation = True
                            print("Found acid to acid chloride activation")
                            break

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_acid_activation
