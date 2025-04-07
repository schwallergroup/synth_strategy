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
    This function detects if the synthetic route involves reduction of an alkyne to an alkane.
    """
    alkyne_reduction_found = False

    def dfs_traverse(node):
        nonlocal alkyne_reduction_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for alkyne in reactants
                alkyne_pattern = Chem.MolFromSmarts("[C]#[C]")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    if reactant_mol.HasSubstructMatch(alkyne_pattern):
                        # Check if the product has the alkyne converted to alkane
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            # If product doesn't have alkyne but has similar structure,
                            # it's likely a reduction occurred
                            if not product_mol.HasSubstructMatch(alkyne_pattern):
                                print("Alkyne reduction detected")
                                alkyne_reduction_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return alkyne_reduction_found
