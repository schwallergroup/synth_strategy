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
    This function detects a linear synthesis strategy that preserves specific functional groups
    (thioether and aryl chloride) throughout the synthesis.
    """
    # Track if we've found the functional groups in the final product
    found_thioether = False
    found_aryl_chloride = False

    # Track if synthesis is linear (no convergent steps)
    is_linear = True

    def dfs_traverse(node):
        nonlocal found_thioether, found_aryl_chloride, is_linear

        if node["type"] == "mol" and not node.get("in_stock", False):
            # Check final product for functional groups
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                thioether_pattern = Chem.MolFromSmarts("[#6]-[#16]-[#6]")
                aryl_chloride_pattern = Chem.MolFromSmarts("c-[Cl]")

                if mol.HasSubstructMatch(thioether_pattern):
                    found_thioether = True

                if mol.HasSubstructMatch(aryl_chloride_pattern):
                    found_aryl_chloride = True

        elif node["type"] == "reaction":
            # Check if reaction is convergent (more than 2 reactants)
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            if len(reactants) > 2:
                is_linear = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if synthesis is linear and preserves both functional groups
    return is_linear and found_thioether and found_aryl_chloride
