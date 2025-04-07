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
    This function detects a late-stage Suzuki coupling of two complex fragments.
    """
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction" and depth <= 1:  # Late stage (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Suzuki coupling pattern
                reactant_mols = [
                    Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
                ]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and len(reactant_mols) >= 2:
                    # Check for boronic acid/ester in one reactant
                    boronic_pattern = Chem.MolFromSmarts("[c]-[B]")
                    # Check for aryl halide in another reactant
                    halide_pattern = Chem.MolFromSmarts("[c]-[Br,I,Cl]")

                    has_boronic = any(
                        r and r.HasSubstructMatch(boronic_pattern)
                        for r in reactant_mols
                    )
                    has_halide = any(
                        r and r.HasSubstructMatch(halide_pattern) for r in reactant_mols
                    )

                    # Check if product has a new biaryl bond
                    if has_boronic and has_halide:
                        # Check if reactants are complex (>15 heavy atoms)
                        complex_fragments = sum(
                            1
                            for r in reactant_mols
                            if r and Descriptors.HeavyAtomCount(r) > 15
                        )

                        if complex_fragments >= 1:
                            print(
                                "Found late-stage Suzuki coupling of complex fragments"
                            )
                            found_pattern = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_pattern
