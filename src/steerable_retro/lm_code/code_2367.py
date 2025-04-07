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
    This function detects if the synthesis includes a tertiary amine formation step.
    """
    tertiary_amine_pattern = Chem.MolFromSmarts("[#7]([#6])([#6])[#6]")
    secondary_amine_pattern = Chem.MolFromSmarts("[#7;H1]([#6])[#6]")
    tertiary_amine_formation_found = False

    def dfs_traverse(node):
        nonlocal tertiary_amine_formation_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            reactant_mol = Chem.MolFromSmiles(reactants)
            product_mol = Chem.MolFromSmiles(product)

            if (
                reactant_mol
                and product_mol
                and reactant_mol.HasSubstructMatch(secondary_amine_pattern)
                and product_mol.HasSubstructMatch(tertiary_amine_pattern)
            ):
                tertiary_amine_formation_found = True
                print("Tertiary amine formation detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return tertiary_amine_formation_found
