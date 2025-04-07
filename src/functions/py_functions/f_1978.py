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
    Detects if the synthesis route involves conversion of a benzyl alcohol to a benzyl chloride
    """
    has_conversion = False

    def dfs_traverse(node):
        nonlocal has_conversion

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for benzyl alcohol to benzyl chloride conversion
                benzyl_alcohol_patt = Chem.MolFromSmarts("[c][CH2][OH]")
                benzyl_chloride_patt = Chem.MolFromSmarts("[c][CH2][Cl]")

                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(benzyl_chloride_patt):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            benzyl_alcohol_patt
                        ):
                            has_conversion = True
                            print("Found alcohol to chloride conversion")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_conversion
