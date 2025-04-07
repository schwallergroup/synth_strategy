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
    This function detects ester reduction to primary alcohol.
    """
    ester_reduction_found = False

    def dfs_traverse(node):
        nonlocal ester_reduction_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if reactant has ester and product has primary alcohol
            reactant_mol = Chem.MolFromSmiles(reactants_part)
            product_mol = Chem.MolFromSmiles(product_part)

            if reactant_mol and product_mol:
                ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]")
                alcohol_pattern = Chem.MolFromSmarts("[CH2][OH]")

                if reactant_mol.HasSubstructMatch(
                    ester_pattern
                ) and product_mol.HasSubstructMatch(alcohol_pattern):
                    print("Found ester reduction to primary alcohol")
                    ester_reduction_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return ester_reduction_found
