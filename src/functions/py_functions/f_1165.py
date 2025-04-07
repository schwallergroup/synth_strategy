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
    This function detects if the synthesis involves thioether formation from an alkyl halide and thiol.
    """
    thioether_formation_detected = False

    def dfs_traverse(node):
        nonlocal thioether_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactants contain alkyl halide and thiol
                alkyl_halide_pattern = Chem.MolFromSmarts("[#6]-[#35,#53,#17]")
                thiol_pattern = Chem.MolFromSmarts("[#16H]")

                # Check if product contains thioether
                thioether_pattern = Chem.MolFromSmarts("[#6]-[#16]-[#6]")

                has_alkyl_halide = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(alkyl_halide_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r) is not None
                )

                has_thiol = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(thiol_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r) is not None
                )

                product_mol = Chem.MolFromSmiles(product)
                has_thioether = (
                    product_mol is not None
                    and product_mol.HasSubstructMatch(thioether_pattern)
                )

                if has_alkyl_halide and has_thiol and has_thioether:
                    print("Detected thioether formation from alkyl halide and thiol")
                    thioether_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return thioether_formation_detected
