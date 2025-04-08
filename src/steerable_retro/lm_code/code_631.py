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
    This function detects a linear synthesis strategy involving nitrile-containing
    intermediates in the synthetic pathway.
    """
    nitrile_intermediates = 0
    linear_synthesis = True

    def dfs_traverse(node):
        nonlocal nitrile_intermediates, linear_synthesis

        if node["type"] == "mol" and not node.get("in_stock", False):
            # Check for nitrile group in intermediates
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[C]#[N]")):
                nitrile_intermediates += 1
                print(f"Nitrile intermediate detected: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check if more than 2 reactants (non-linear synthesis)
            if len(reactants) > 2:
                linear_synthesis = False
                print(f"Non-linear step detected with {len(reactants)} reactants")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if it's a linear synthesis with at least one nitrile intermediate
    return linear_synthesis and nitrile_intermediates > 0
