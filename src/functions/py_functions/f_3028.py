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
    This function detects if the synthetic route uses late-stage amide formation.
    """
    amide_formation_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for acyl chloride and aniline in reactants
            acyl_chloride_pattern = Chem.MolFromSmarts("[Cl]-[C](=O)-[c]")
            aniline_pattern = Chem.MolFromSmarts("[NH2]-[c]")
            amide_pattern = Chem.MolFromSmarts("[NH]-[C](=O)-[c]")

            has_acyl_chloride = False
            has_aniline = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    if reactant_mol.HasSubstructMatch(acyl_chloride_pattern):
                        has_acyl_chloride = True
                    if reactant_mol.HasSubstructMatch(aniline_pattern):
                        has_aniline = True

            product_mol = Chem.MolFromSmiles(product)
            if (
                has_acyl_chloride
                and has_aniline
                and product_mol
                and product_mol.HasSubstructMatch(amide_pattern)
            ):
                amide_formation_depth = depth
                print(f"Amide formation detected at depth {depth}: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it late-stage if it's in the first half of the synthesis (lower depth in retrosynthesis)
    is_late_stage = (
        amide_formation_depth is not None and amide_formation_depth <= max_depth / 2
    )

    print(
        f"Amide formation depth: {amide_formation_depth}, Max depth: {max_depth}, Is late-stage: {is_late_stage}"
    )
    return is_late_stage
