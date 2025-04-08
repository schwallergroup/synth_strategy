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
    Detects if the synthesis route includes a late-stage nitro reduction to amine.
    Late stage is defined as occurring at depth 0 or 1.
    """
    nitro_reduction_found = False
    max_depth_for_late_stage = 1

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_found

        if node["type"] == "reaction" and depth <= max_depth_for_late_stage:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a nitro reduction reaction
            reactant_mol = Chem.MolFromSmiles(reactants[0])
            product_mol = Chem.MolFromSmiles(product)

            if reactant_mol and product_mol:
                nitro_pattern = Chem.MolFromSmarts("[c][N+](=O)[O-]")
                amine_pattern = Chem.MolFromSmarts("[c][NH2]")

                if reactant_mol.HasSubstructMatch(nitro_pattern) and product_mol.HasSubstructMatch(
                    amine_pattern
                ):
                    # Confirm nitro group is reduced to amine
                    nitro_matches = reactant_mol.GetSubstructMatches(nitro_pattern)
                    for match in nitro_matches:
                        aromatic_carbon_idx = match[0]
                        if not product_mol.HasSubstructMatch(nitro_pattern) or len(
                            product_mol.GetSubstructMatches(nitro_pattern)
                        ) < len(nitro_matches):
                            nitro_reduction_found = True
                            print(f"Late-stage nitro reduction detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return nitro_reduction_found
