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
    This function detects if a synthetic route contains a late-stage fragment coupling
    via Williamson ether synthesis (in the first half of the synthesis).
    """
    found_williamson = False
    max_depth = 0
    williamson_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal found_williamson, max_depth, williamson_depth

        # Track maximum depth to determine what's "late stage"
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Williamson ether synthesis
            # Look for phenol/alcohol + benzyl bromide → ether
            phenol_pattern = Chem.MolFromSmarts("[c,C][OH]")
            bromide_pattern = Chem.MolFromSmarts("[c,C][C][Br]")
            ether_pattern = Chem.MolFromSmarts("[c,C][O][C]")

            has_phenol = False
            has_bromide = False

            for reactant in reactants:
                if reactant.strip():
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(phenol_pattern):
                                has_phenol = True
                            if mol.HasSubstructMatch(bromide_pattern):
                                has_bromide = True
                    except:
                        continue

            has_ether = False
            if product.strip():
                try:
                    mol = Chem.MolFromSmiles(product)
                    if mol and mol.HasSubstructMatch(ether_pattern):
                        has_ether = True
                except:
                    pass

            # If transformation matches Williamson pattern
            if has_phenol and has_bromide and has_ether:
                found_williamson = True
                williamson_depth = depth
                print(f"Found Williamson ether synthesis at depth {depth}: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if Williamson ether synthesis is in the first half (late stage)
    is_late_stage = williamson_depth is not None and williamson_depth <= max_depth / 2

    result = found_williamson and is_late_stage
    print(
        f"Late-stage fragment coupling strategy detected: {result} (depth: {williamson_depth}, max depth: {max_depth})"
    )
    return result
