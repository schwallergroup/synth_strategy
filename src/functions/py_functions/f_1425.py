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
    This function detects a synthetic strategy involving C-O ether bond formation
    between an alkyl bromide and bromophenol.
    """
    found_ether_formation = False

    def dfs_traverse(node):
        nonlocal found_ether_formation

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for bromophenol
            bromophenol_pattern = Chem.MolFromSmarts("[OH][c]1[c]([Br])[c][c][c][c]1")

            # Check for alkyl bromide
            alkyl_bromide_pattern = Chem.MolFromSmarts("[C][Br]")

            # Check for ether in product
            ether_pattern = Chem.MolFromSmarts("[C][O][c]1[c]([Br])[c][c][c][c]1")

            reactant_has_bromophenol = False
            reactant_has_alkyl_bromide = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(bromophenol_pattern):
                        reactant_has_bromophenol = True
                    if mol.HasSubstructMatch(alkyl_bromide_pattern):
                        reactant_has_alkyl_bromide = True

            product_mol = Chem.MolFromSmiles(product)
            if (
                reactant_has_bromophenol
                and reactant_has_alkyl_bromide
                and product_mol
                and product_mol.HasSubstructMatch(ether_pattern)
            ):
                found_ether_formation = True
                print("Found ether formation between alkyl bromide and bromophenol")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    if found_ether_formation:
        print("Detected ether formation with bromophenol strategy")

    return found_ether_formation
