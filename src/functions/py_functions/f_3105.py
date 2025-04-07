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
    Detects a synthetic strategy involving the conversion of an alcohol to a bromide
    as an activation step for subsequent alkylation.
    """
    has_alcohol_to_bromide = False

    def dfs_traverse(node):
        nonlocal has_alcohol_to_bromide

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alcohol to bromide conversion
            # Pattern: [C][OH] → [C][Br]

            for reactant in reactants:
                r_mol = Chem.MolFromSmiles(reactant)
                if not r_mol:
                    continue

                # Check for alcohol pattern
                alcohol_patt = Chem.MolFromSmarts("[C][OH]")
                if r_mol.HasSubstructMatch(alcohol_patt):
                    # Check if product has C-Br
                    p_mol = Chem.MolFromSmiles(product)
                    if not p_mol:
                        continue

                    bromide_patt = Chem.MolFromSmarts("[C][Br]")
                    if p_mol.HasSubstructMatch(bromide_patt):
                        has_alcohol_to_bromide = True
                        print("Detected alcohol to bromide conversion")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return has_alcohol_to_bromide
