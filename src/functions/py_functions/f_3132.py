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
    This function detects if the synthesis introduces a nitroimidazole group in the late stage.
    """
    nitroimidazole_introduced_late = False

    def dfs_traverse(node, depth=0):
        nonlocal nitroimidazole_introduced_late

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Nitroimidazole pattern
                nitroimidazole_pattern = Chem.MolFromSmarts("[n]1cnc([N+](=O)[O-])c1")

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(
                        nitroimidazole_pattern
                    ):
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                nitroimidazole_pattern
                            ):
                                nitroimidazole_introduced_late = True
                                print(
                                    f"Nitroimidazole introduction detected at depth {depth}"
                                )
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return nitroimidazole_introduced_late
