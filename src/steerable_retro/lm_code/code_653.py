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
    This function detects if a chloro-aromatic fragment is introduced in the synthesis.
    """
    chloro_aromatic_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal chloro_aromatic_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for chloro-aromatic in reactants
                chloro_aromatic_pattern = Chem.MolFromSmarts("[c][Cl]")

                for reactant in reactants:
                    if (
                        reactant
                        and Chem.MolFromSmiles(reactant)
                        and Chem.MolFromSmiles(reactant).HasSubstructMatch(chloro_aromatic_pattern)
                    ):
                        # Check if the product also has the chloro-aromatic
                        if (
                            product
                            and Chem.MolFromSmiles(product)
                            and Chem.MolFromSmiles(product).HasSubstructMatch(
                                chloro_aromatic_pattern
                            )
                        ):
                            print("Chloro-aromatic fragment introduction detected at depth", depth)
                            chloro_aromatic_detected = True
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    return chloro_aromatic_detected
