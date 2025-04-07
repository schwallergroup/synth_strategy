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
    This function detects a strategy where two heteroaromatic fragments are coupled via sulfonamide formation.
    """
    has_sulfonamide_formation = False

    def dfs_traverse(node):
        nonlocal has_sulfonamide_formation

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for sulfonamide formation between heteroaromatic fragments
            has_sulfonyl_chloride = any(
                Chem.MolFromSmiles(r)
                and Chem.MolFromSmiles(r).HasSubstructMatch(
                    Chem.MolFromSmarts("[S](=[O])(=[O])[Cl]")
                )
                for r in reactants
            )
            has_amine = any(
                Chem.MolFromSmiles(r)
                and Chem.MolFromSmiles(r).HasSubstructMatch(
                    Chem.MolFromSmarts("[NH2][c]")
                )
                for r in reactants
            )
            has_sulfonamide_product = Chem.MolFromSmiles(
                product
            ) and Chem.MolFromSmiles(product).HasSubstructMatch(
                Chem.MolFromSmarts("[NH][S](=[O])(=[O])[c]")
            )

            if has_sulfonyl_chloride and has_amine and has_sulfonamide_product:
                has_sulfonamide_formation = True
                print("Found sulfonamide formation between heteroaromatic fragments")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return has_sulfonamide_formation
