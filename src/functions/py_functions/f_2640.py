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
    Detects if the synthesis involves conversion of an alcohol to a chloride.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                if len(reactants) >= 1:
                    product_mol = Chem.MolFromSmiles(product)

                    # Check for benzyl chloride pattern in product
                    benzyl_chloride_pattern = Chem.MolFromSmarts("[c][C][Cl]")

                    if product_mol and product_mol.HasSubstructMatch(
                        benzyl_chloride_pattern
                    ):
                        # Check if reactant has benzyl alcohol
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            benzyl_alcohol_pattern = Chem.MolFromSmarts("[c][C][OH]")

                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                benzyl_alcohol_pattern
                            ):
                                print("Detected alcohol to chloride conversion")
                                result = True
                                break
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return result
