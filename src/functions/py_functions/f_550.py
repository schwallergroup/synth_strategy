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
    This function detects if the synthesis includes formamide formation from amines.
    """
    formamide_formations = 0

    def dfs_traverse(node):
        nonlocal formamide_formations

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amine to formamide conversion
            has_amine = False
            for reactant in reactants:
                if reactant:
                    r_mol = Chem.MolFromSmiles(reactant)
                    if r_mol and r_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[NH2][c]")
                    ):
                        has_amine = True

            if has_amine and product:
                p_mol = Chem.MolFromSmiles(product)
                if p_mol and p_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH][CH]=O")):
                    print("Formamide formation detected")
                    formamide_formations += 1

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return formamide_formations >= 1
