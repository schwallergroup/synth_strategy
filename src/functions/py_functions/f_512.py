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
    This function detects late-stage amide bond formation.
    """
    amide_formation_depth = None
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide formation
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
                if product_mol.HasSubstructMatch(amide_pattern):
                    # Check if any reactant doesn't have the amide
                    amide_in_all_reactants = True
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if not reactant_mol.HasSubstructMatch(amide_pattern):
                                amide_in_all_reactants = False
                                break

                    if not amide_in_all_reactants:
                        # This reaction forms an amide
                        if (
                            amide_formation_depth is None
                            or depth < amide_formation_depth
                        ):
                            amide_formation_depth = depth
                            print(
                                f"Amide formation detected at depth {depth} in reaction: {rsmi}"
                            )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it late-stage if it's in the first half of the synthesis
    # (remember depth 0 is the final step)
    if amide_formation_depth is not None and amide_formation_depth <= max_depth / 2:
        print(
            f"Late-stage amide formation detected at depth {amide_formation_depth} (max depth: {max_depth})"
        )
        return True
    return False
