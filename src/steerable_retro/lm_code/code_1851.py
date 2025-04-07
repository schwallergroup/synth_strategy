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
    Detects a strategy involving early-stage biaryl C-C bond formation.
    """
    biaryl_coupling_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal biaryl_coupling_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for biaryl coupling - tertiary carbon connected to aromatic ring
                if len(reactants) >= 2:  # Need at least two fragments
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                    product_mol = Chem.MolFromSmiles(product)

                    # Patterns for tertiary carbon and aromatic ring
                    tert_carbon_pattern = Chem.MolFromSmarts("[C]([C])([C])[C]")
                    aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")

                    # Check if reactants have these patterns separately and product has them connected
                    has_tert_carbon = any(
                        r and r.HasSubstructMatch(tert_carbon_pattern) for r in reactant_mols if r
                    )
                    has_aromatic = any(
                        r and r.HasSubstructMatch(aromatic_pattern) for r in reactant_mols if r
                    )

                    if has_tert_carbon and has_aromatic:
                        # Check if product has tertiary carbon connected to aromatic ring
                        biaryl_pattern = Chem.MolFromSmarts("[C]([C])([C])[c]")
                        if product_mol and product_mol.HasSubstructMatch(biaryl_pattern):
                            biaryl_coupling_depth = depth
                            print(f"Found biaryl coupling at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if biaryl coupling occurred in early stage (second half of synthesis)
    if biaryl_coupling_depth is not None and biaryl_coupling_depth > max_depth / 2:
        print(
            f"Confirmed early-stage biaryl coupling at depth {biaryl_coupling_depth} (max depth: {max_depth})"
        )
        return True

    return False
