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
    This function detects if the synthesis includes a benzylic alcohol oxidation step.
    """
    found_benzylic_oxidation = False

    def dfs_traverse(node):
        nonlocal found_benzylic_oxidation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for benzylic alcohol oxidation
                benzylic_alcohol_pattern = Chem.MolFromSmarts("[c][C][OH]")
                methylene_pattern = Chem.MolFromSmarts("[c][C][c]")

                for reactant in reactants:
                    try:
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol and r_mol.HasSubstructMatch(benzylic_alcohol_pattern):
                            p_mol = Chem.MolFromSmiles(product)
                            if p_mol and p_mol.HasSubstructMatch(methylene_pattern):
                                found_benzylic_oxidation = True
                                print("Found benzylic alcohol oxidation")
                    except:
                        continue

        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Benzylic alcohol oxidation strategy detected: {found_benzylic_oxidation}")
    return found_benzylic_oxidation
