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
    This function detects if the synthetic route involves phenol protection with a carbonate group.
    """
    phenol_protection_found = False

    def dfs_traverse(node):
        nonlocal phenol_protection_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a phenol protection reaction
                try:
                    # Look for phenol in reactants
                    phenol_pattern = Chem.MolFromSmarts("[OH][c]")
                    carbonate_pattern = Chem.MolFromSmarts("[O][C](=[O])[O][C]")

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            phenol_pattern
                        ):
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol and product_mol.HasSubstructMatch(
                                carbonate_pattern
                            ):
                                print(f"Found phenol protection: {rsmi}")
                                phenol_protection_found = True
                except Exception as e:
                    print(f"Error in phenol protection detection: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return phenol_protection_found
