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
    Detects if the synthesis uses a late-stage reductive amination strategy
    (in the final 2-3 steps of synthesis)
    """
    reductive_amination_found = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal reductive_amination_found, max_depth

        max_depth = max(max_depth, depth)

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a reductive amination
            # Look for aldehyde pattern in reactants
            aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
            amine_pattern = Chem.MolFromSmarts("[NX3;!$(NC=O)]")

            has_aldehyde = False
            has_amine = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(aldehyde_pattern):
                        has_aldehyde = True
                    if mol.HasSubstructMatch(amine_pattern):
                        has_amine = True

            # Check if product has a new C-N bond
            if has_aldehyde and has_amine and depth <= 2:  # Late stage (depth 0-2)
                reductive_amination_found = True
                print(f"Found late-stage reductive amination at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return reductive_amination_found
