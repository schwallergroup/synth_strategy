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
    This function detects if the synthetic route contains multiple sulfonamide formation reactions.
    """
    sulfonamide_formation_count = 0

    def dfs_traverse(node):
        nonlocal sulfonamide_formation_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Check if this is a sulfonamide formation reaction
                # Look for sulfonyl chloride in reactants
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants.split(".")]
                product_mol = Chem.MolFromSmiles(product)

                sulfonyl_chloride_pattern = Chem.MolFromSmarts("[#16](=[O])(=[O])[Cl]")
                amine_pattern = Chem.MolFromSmarts("[N;!$(N=*);!$(N#*)]")
                sulfonamide_pattern = Chem.MolFromSmarts("[N][S](=[O])(=[O])")

                has_sulfonyl_chloride = any(
                    mol.HasSubstructMatch(sulfonyl_chloride_pattern)
                    for mol in reactant_mols
                    if mol
                )
                has_amine = any(
                    mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols if mol
                )
                has_sulfonamide_product = product_mol and product_mol.HasSubstructMatch(
                    sulfonamide_pattern
                )

                if has_sulfonyl_chloride and has_amine and has_sulfonamide_product:
                    sulfonamide_formation_count += 1
                    print(f"Found sulfonamide formation reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    result = sulfonamide_formation_count >= 2
    print(
        f"Multiple sulfonamide formations detected: {result} (count: {sulfonamide_formation_count})"
    )
    return result
