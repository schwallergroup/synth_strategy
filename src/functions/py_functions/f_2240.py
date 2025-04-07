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
    This function detects if the synthesis involves formation of a urea linker
    connecting two aromatic systems.
    """
    urea_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal urea_formation_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for urea pattern in product
            urea_pattern = Chem.MolFromSmarts("[c,n][NH][C](=[O])[NH][c,n]")
            product_mol = Chem.MolFromSmiles(product)

            if product_mol and product_mol.HasSubstructMatch(urea_pattern):
                # Check if urea wasn't present in reactants
                urea_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(urea_pattern):
                        urea_in_reactants = True
                        break

                if not urea_in_reactants:
                    urea_formation_detected = True
                    print(f"Urea linker formation detected at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return urea_formation_detected
