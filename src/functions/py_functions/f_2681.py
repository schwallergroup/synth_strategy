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
    This function detects ketone reduction to alcohol in the synthesis.
    """
    ketone_reduction_found = False

    def dfs_traverse(node):
        nonlocal ketone_reduction_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ketone reduction pattern
            if len(reactants) > 0:
                reactant_mol = Chem.MolFromSmiles(reactants[0])
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    # Check for ketone in reactant
                    ketone_pattern = Chem.MolFromSmarts("[C]=O")
                    # Check for alcohol in product
                    alcohol_pattern = Chem.MolFromSmarts("[C][OH]")

                    if (
                        reactant_mol.HasSubstructMatch(ketone_pattern)
                        and product_mol.HasSubstructMatch(alcohol_pattern)
                        and not reactant_mol.HasSubstructMatch(alcohol_pattern)
                    ):
                        ketone_reduction_found = True
                        print("Ketone reduction detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return ketone_reduction_found
