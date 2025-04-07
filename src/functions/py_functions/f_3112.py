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
    This function detects early-stage heterocycle formation,
    specifically looking for benzofuranone formation.
    """
    benzofuranone_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c]2[c]1[n][c](=O)[o]2")
    ring_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal ring_formation_detected

        if node["type"] == "reaction" and depth >= 3:  # Early stage (high depth)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains benzofuranone but reactants don't
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(benzofuranone_pattern):
                    reactant_has_pattern = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            benzofuranone_pattern
                        ):
                            reactant_has_pattern = True
                            break

                    if not reactant_has_pattern:
                        print(
                            f"Early-stage heterocycle formation detected at depth {depth}"
                        )
                        ring_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return ring_formation_detected
