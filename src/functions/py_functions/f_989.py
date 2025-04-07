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
    This function detects if the synthetic route contains O-alkylation of a phenol.
    """
    found_o_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_o_alkylation

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for phenol in reactants
                phenol_pattern = Chem.MolFromSmarts("c[OH]")

                # Check for alkoxy group in product
                alkoxy_pattern = Chem.MolFromSmarts("c[O][C]")

                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(alkoxy_pattern):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            phenol_pattern
                        ):
                            print(f"Found O-alkylation of phenol at depth {depth}")
                            found_o_alkylation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return found_o_alkylation
