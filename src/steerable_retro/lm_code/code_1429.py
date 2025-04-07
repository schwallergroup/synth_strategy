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
    This function detects a synthetic strategy involving late-stage
    heterocycle (thiazole) formation.
    """
    # Track depth of thiazole formation
    thiazole_formation_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal thiazole_formation_depth, max_depth

        # Update max depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for thiazole formation
                thiazole_pattern = Chem.MolFromSmarts("c1scnc1")

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(thiazole_pattern):
                        # Check if any reactant doesn't have thiazole
                        has_thiazole_in_reactants = False
                        for reactant in reactants:
                            try:
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol and reactant_mol.HasSubstructMatch(
                                    thiazole_pattern
                                ):
                                    has_thiazole_in_reactants = True
                                    break
                            except:
                                continue

                        if not has_thiazole_in_reactants:
                            thiazole_formation_depth = depth
                            print(f"Found thiazole formation at depth {depth}")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it late-stage if it occurs in the first half of the synthesis
    # (remember: lower depth = later in synthesis)
    if thiazole_formation_depth is not None and thiazole_formation_depth <= max_depth / 2:
        print(
            f"Confirmed late-stage heterocycle formation (depth {thiazole_formation_depth} out of max {max_depth})"
        )
        return True
    return False
