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
    This function detects a Suzuki coupling reaction in the synthetic route.
    """
    found_suzuki = False

    # SMARTS patterns for Suzuki coupling components
    boronic_pattern = Chem.MolFromSmarts("[#6]~[#5](~[#8])~[#8]")  # Boronic acid/ester
    aryl_halide_pattern = Chem.MolFromSmarts("[#6]~[#53,#35,#17,#9]")  # Aryl halide

    def dfs_traverse(node, depth=0):
        nonlocal found_suzuki

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                reactants = [
                    Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r
                ]
                product = Chem.MolFromSmiles(product_str) if product_str else None

                if product and len(reactants) >= 2:
                    # Check for boronic acid/ester and aryl halide in reactants
                    has_boronic = any(
                        r and r.HasSubstructMatch(boronic_pattern) for r in reactants
                    )
                    has_aryl_halide = any(
                        r and r.HasSubstructMatch(aryl_halide_pattern)
                        for r in reactants
                    )

                    if has_boronic and has_aryl_halide:
                        found_suzuki = True
                        print(f"Found Suzuki coupling at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Suzuki coupling strategy detected: {found_suzuki}")
    return found_suzuki
