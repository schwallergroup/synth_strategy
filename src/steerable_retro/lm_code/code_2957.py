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
    This function detects if the synthetic route contains a late-stage sulfonylation (depth 0 or 1).
    """
    late_stage_sulfonylation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_sulfonylation_found

        if node["type"] == "reaction" and depth <= 1:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants.split(".")]
                product_mol = Chem.MolFromSmiles(product)

                # Patterns for sulfonylation
                sulfonyl_chloride_pattern = Chem.MolFromSmarts("[#16](=[O])(=[O])[Cl]")
                sulfonamide_pattern = Chem.MolFromSmarts("[N][S](=[O])(=[O])")

                has_sulfonyl_chloride = any(
                    mol and mol.HasSubstructMatch(sulfonyl_chloride_pattern)
                    for mol in reactant_mols
                    if mol
                )
                has_sulfonamide_product = product_mol and product_mol.HasSubstructMatch(
                    sulfonamide_pattern
                )

                if has_sulfonyl_chloride and has_sulfonamide_product:
                    late_stage_sulfonylation_found = True
                    print(f"Found late-stage sulfonylation at depth {depth}: {rsmi}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage sulfonylation detected: {late_stage_sulfonylation_found}")
    return late_stage_sulfonylation_found
