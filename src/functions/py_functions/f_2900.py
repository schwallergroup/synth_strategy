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
    Detects the use of Stille coupling (organotin reagent) for C-C bond formation
    """
    found_stille_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_stille_coupling

        if node["type"] == "reaction":
            # Extract reaction information
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for organotin reagent in reactants
            organotin_pattern = r"Sn"
            vinyl_pattern = Chem.MolFromSmarts("[#6]=[#6]")

            # Look for organotin in reactants and vinyl in product
            has_organotin = any(organotin_pattern in r for r in reactants)

            if has_organotin:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(vinyl_pattern):
                    found_stille_coupling = True
                    print(f"Found Stille coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Stille coupling strategy detected: {found_stille_coupling}")
    return found_stille_coupling
