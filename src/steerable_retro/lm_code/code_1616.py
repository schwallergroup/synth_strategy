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
    This function detects if the route contains a late-stage fluorine displacement,
    typically to form an ether.
    """
    found = False

    def dfs_traverse(node):
        nonlocal found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for fluorine in reactants
            reactant_has_fluorine = False
            for reactant in reactants:
                if reactant:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]-[#9]")):
                        reactant_has_fluorine = True

            # Check for ether in product
            product_mol = Chem.MolFromSmiles(product)
            product_has_ether = product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts("[#6]-[#8]-[#6]")
            )

            if reactant_has_fluorine and product_has_ether:
                # Check if this is a late-stage reaction (depth <= 1)
                depth = node.get("metadata", {}).get("depth", -1)
                if depth <= 1:  # Assuming depth 0 or 1 is late-stage
                    found = True
                    print(f"Found late-stage fluorine displacement: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage fluorine displacement: {found}")
    return found
