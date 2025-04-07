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
    This function detects a synthetic strategy using a bifunctional linker (like dibromoethane)
    for fragment connection.
    """
    bifunctional_linker_found = False

    def dfs_traverse(node):
        nonlocal bifunctional_linker_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for dibromoethane or similar bifunctional linkers
            dibromo_pattern = Chem.MolFromSmarts("[Br,Cl,I][#6][#6][Br,Cl,I]")

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if (
                    reactant_mol
                    and dibromo_pattern
                    and reactant_mol.HasSubstructMatch(dibromo_pattern)
                ):
                    # Check if product has incorporated one end of the linker
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[#8,#7,#16][#6][#6][Br,Cl,I]")
                    ):
                        print("Found bifunctional linker (dibromoethane) incorporation")
                        bifunctional_linker_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return bifunctional_linker_found
