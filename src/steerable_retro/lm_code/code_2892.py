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
    This function detects sulfonamide formation from an amine and a sulfonyl chloride.
    """
    sulfonamide_formed = False

    def dfs_traverse(node):
        nonlocal sulfonamide_formed

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for sulfonamide formation
            sulfonamide_pattern = Chem.MolFromSmarts("[S](=[O])(=[O])[N]")
            amine_pattern = Chem.MolFromSmarts("[NH2]")
            sulfonyl_chloride_pattern = Chem.MolFromSmarts("[S](=[O])(=[O])[Cl]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and product_mol.HasSubstructMatch(sulfonamide_pattern):
                has_amine = any(m and m.HasSubstructMatch(amine_pattern) for m in reactant_mols)
                has_sulfonyl_chloride = any(
                    m and m.HasSubstructMatch(sulfonyl_chloride_pattern) for m in reactant_mols
                )

                if has_amine and has_sulfonyl_chloride:
                    sulfonamide_formed = True
                    print("Found sulfonamide formation from amine and sulfonyl chloride")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return sulfonamide_formed
