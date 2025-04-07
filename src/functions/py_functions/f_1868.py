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
    This function detects if the synthesis involves late-stage amide formation.
    """
    amide_formation_depth = None
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_depth, max_depth

        max_depth = max(max_depth, depth)

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Amide pattern
            amide_pattern = Chem.MolFromSmarts("[N;X3]([C;X3]=[O;X1])")

            product_mol = Chem.MolFromSmiles(product)

            if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                # Check if amide was not present in reactants
                amide_in_reactants = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(amide_pattern):
                        amide_in_reactants = True
                        break

                if not amide_in_reactants:
                    amide_formation_depth = depth
                    print(f"Amide formation detected at depth {depth}: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Consider it late-stage if it's in the first third of the synthesis
    if amide_formation_depth is not None and amide_formation_depth <= max_depth // 3:
        print(
            f"Late-stage amide formation confirmed at depth {amide_formation_depth} (max depth: {max_depth})"
        )
        return True
    return False
