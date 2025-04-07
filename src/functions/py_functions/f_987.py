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
    This function detects if the synthetic route contains a late-stage nitro reduction
    (nitro group reduced to amine in the final or penultimate step).
    """
    found_nitro_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_nitro_reduction

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Final or penultimate step (depth 0 or 1)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Check if reactant contains nitro group and product contains amine
                reactant_mol = Chem.MolFromSmiles(reactants)
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
                    amine_pattern = Chem.MolFromSmarts("[NH2]")

                    if reactant_mol.HasSubstructMatch(
                        nitro_pattern
                    ) and product_mol.HasSubstructMatch(amine_pattern):
                        print(f"Found nitro reduction at depth {depth}")
                        found_nitro_reduction = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return found_nitro_reduction
