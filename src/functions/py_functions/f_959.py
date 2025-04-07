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
    This function detects the use of silyl protection in a synthetic route.
    Specifically, it looks for the formation of O-Si bonds in protection steps.
    """
    silyl_protection_found = False

    def dfs_traverse(node):
        nonlocal silyl_protection_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for silyl protection: formation of O-Si bond
            # Look for patterns like [O:n][Si:m] in the product but not in reactants
            o_si_pattern = Chem.MolFromSmarts("[#8]-[#14]")

            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(o_si_pattern):
                # Check if this bond was formed in this reaction
                reactant_has_o_si = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(o_si_pattern):
                        reactant_has_o_si = True
                        break

                if not reactant_has_o_si:
                    print("Silyl protection detected: O-Si bond formed")
                    silyl_protection_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return silyl_protection_found
