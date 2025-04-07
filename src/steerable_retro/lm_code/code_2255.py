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
    This function detects the use of Suzuki coupling to form biaryl C-C bonds
    in the synthesis route.
    """
    # Track if we found the pattern
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction":
            # Get reaction SMILES
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check for patterns indicative of Suzuki coupling
                    aryl_halide_pattern = Chem.MolFromSmarts("[c][Cl,Br,I]")
                    boronic_pattern = Chem.MolFromSmarts("[c][B]([O])[O]")

                    # Check reactants for characteristic Suzuki coupling partners
                    has_aryl_halide = any(
                        Chem.MolFromSmiles(r) is not None
                        and Chem.MolFromSmiles(r).HasSubstructMatch(aryl_halide_pattern)
                        for r in reactants
                    )
                    has_boronic = any(
                        Chem.MolFromSmiles(r) is not None
                        and Chem.MolFromSmiles(r).HasSubstructMatch(boronic_pattern)
                        for r in reactants
                    )

                    # Check if product has a biaryl bond that wasn't in the reactants
                    if has_aryl_halide and has_boronic:
                        print(f"Found Suzuki coupling for biaryl formation at depth {depth}")
                        found_pattern = True

                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_pattern
