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
    This function detects if the synthetic route involves a convergent approach
    with fragment coupling via C-N bond formation.
    """
    # Track if we found the pattern
    found_fragment_coupling = False

    def dfs_traverse(node):
        nonlocal found_fragment_coupling

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if we have multiple reactants (fragments)
            if len(reactants) >= 2:
                # Check if one fragment has a nitrogen heterocycle
                nitrogen_heterocycle_pattern = "[#7;R]"

                fragment_count = 0
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(
                        Chem.MolFromSmarts(nitrogen_heterocycle_pattern)
                    ):
                        fragment_count += 1

                # If we have at least two fragments with nitrogen heterocycles
                if fragment_count >= 2:
                    # Check if product has both fragments connected
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # Count the number of fragments in the product
                        frags = Chem.GetMolFrags(product_mol)
                        if len(frags) == 1:  # Single connected product
                            print("Found convergent fragment coupling via C-N bond formation")
                            found_fragment_coupling = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_fragment_coupling
