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
    This function detects if the synthetic route contains a late-stage amide formation
    (in the final or penultimate step).
    """
    amide_formation_depths = []
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                reactants = Chem.MolFromSmiles(reactants_smiles)
                product = Chem.MolFromSmiles(product_smiles)

                if reactants and product:
                    # Carboxylic acid or derivative pattern
                    acid_pattern = Chem.MolFromSmarts("[C](=[O])[O,N]")
                    # Amine pattern
                    amine_pattern = Chem.MolFromSmarts("[N;!$(N=*);!$(NC=O)]")
                    # Amide pattern
                    amide_pattern = Chem.MolFromSmarts("[N;!$(N=*)][C](=[O])")

                    # Check for amide formation
                    if (
                        reactants.HasSubstructMatch(acid_pattern)
                        and reactants.HasSubstructMatch(amine_pattern)
                        and product.HasSubstructMatch(amide_pattern)
                    ):
                        amide_formation_depths.append(depth)
                        print(f"Amide formation detected at depth {depth}: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if there's an amide formation in the final or penultimate step
    for depth in amide_formation_depths:
        if depth <= 1:  # Final or penultimate step (depth 0 or 1)
            return True

    return False
