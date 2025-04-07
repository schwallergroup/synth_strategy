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
    This function detects early-stage aromatic hydroxylation.
    """
    hydroxylation_depth = None
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal hydroxylation_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for hydroxylation
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                phenol_pattern = Chem.MolFromSmarts("[c][OH]")
                if product_mol.HasSubstructMatch(phenol_pattern):
                    # Check if any reactant doesn't have the hydroxyl
                    hydroxyl_in_all_reactants = True
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if not reactant_mol.HasSubstructMatch(phenol_pattern):
                                hydroxyl_in_all_reactants = False
                                break

                    if not hydroxyl_in_all_reactants:
                        # This reaction forms a phenol
                        if hydroxylation_depth is None or depth > hydroxylation_depth:
                            hydroxylation_depth = depth
                            print(
                                f"Hydroxylation detected at depth {depth} in reaction: {rsmi}"
                            )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it early-stage if it's in the second half of the synthesis
    # (remember higher depth is earlier in the synthesis)
    if hydroxylation_depth is not None and hydroxylation_depth > max_depth / 2:
        print(
            f"Early-stage hydroxylation detected at depth {hydroxylation_depth} (max depth: {max_depth})"
        )
        return True
    return False
