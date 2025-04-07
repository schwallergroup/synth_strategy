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
    Detects if the synthesis involves conversion of a terminal halide to an azide group.
    """
    has_azide_formation = False

    def dfs_traverse(node):
        nonlocal has_azide_formation

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains azide group
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts("[#6][N]=[N+]=[N-]")
            ):
                # Check if reactants contain alkyl halide
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[#6][Br,Cl,I]")
                    ):
                        has_azide_formation = True
                        print("Detected terminal halide to azide conversion")
                        break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_azide_formation
