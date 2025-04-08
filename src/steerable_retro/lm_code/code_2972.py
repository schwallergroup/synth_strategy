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
    Detects if the route employs a halogen-metal exchange strategy
    (typically bromide to boronic ester) for cross-coupling preparation.
    """
    has_halogen_metal_exchange = False

    def dfs_traverse(node, depth=0):
        nonlocal has_halogen_metal_exchange

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl halide to boronic ester transformation
            aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Br,I]")
            boronic_pattern = Chem.MolFromSmarts("[c]-[B](-[O])-[O]")

            try:
                reactant_mol = Chem.MolFromSmiles(reactants[0])
                product_mol = Chem.MolFromSmiles(product)

                if (
                    reactant_mol
                    and product_mol
                    and reactant_mol.HasSubstructMatch(aryl_halide_pattern)
                    and product_mol.HasSubstructMatch(boronic_pattern)
                ):
                    has_halogen_metal_exchange = True
                    print(f"Found halogen-metal exchange at depth {depth}")
            except:
                pass

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return has_halogen_metal_exchange
