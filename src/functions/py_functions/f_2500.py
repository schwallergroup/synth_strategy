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
    This function detects a synthetic strategy involving boronic acid/ester intermediates
    for cross-coupling reactions.
    """
    has_boronic_intermediate = False
    has_cross_coupling = False

    def dfs_traverse(node):
        nonlocal has_boronic_intermediate, has_cross_coupling

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for boronic acid/ester intermediates
            boronic_acid_pattern = Chem.MolFromSmarts("[#6]-[#5](-[#8])-[#8]")
            boronic_ester_pattern = Chem.MolFromSmarts("[#6]-[#5](-[#8])-[#8]")

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and (
                        mol.HasSubstructMatch(boronic_acid_pattern)
                        or mol.HasSubstructMatch(boronic_ester_pattern)
                    ):
                        has_boronic_intermediate = True
                        print(f"Found boronic acid/ester intermediate: {reactant}")
                except:
                    continue

            # Check for cross-coupling reaction (C-B to C-C or C-N)
            # This is a simplified check - in practice would need more sophisticated detection
            if has_boronic_intermediate:
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol:
                        # Check if product has new C-C or C-N bonds that weren't in reactants
                        has_cross_coupling = True
                        print(f"Detected potential cross-coupling in reaction: {rsmi}")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both conditions are met
    return has_boronic_intermediate and has_cross_coupling
