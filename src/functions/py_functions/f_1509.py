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
    Detects if the synthetic route proceeds through a nitrile intermediate
    that is later transformed into another functional group.
    """
    nitrile_depths = []
    nitrile_transformed = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_depths, nitrile_transformed

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains nitrile
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                nitrile_pattern = Chem.MolFromSmarts("[#6]#[#7]")
                if product_mol.HasSubstructMatch(nitrile_pattern):
                    nitrile_depths.append(depth)
                    print(f"Nitrile formed at depth {depth}")

            # Check if reactants contain nitrile but product doesn't (transformation)
            reactant_has_nitrile = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[#6]#[#7]")
                ):
                    reactant_has_nitrile = True
                    break

            if (
                reactant_has_nitrile
                and product_mol
                and not product_mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]#[#7]"))
            ):
                nitrile_transformed = True
                print(f"Nitrile transformed at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if nitrile was formed and then transformed
    if len(nitrile_depths) > 0 and nitrile_transformed:
        print("Nitrile intermediate strategy detected")
        return True
    return False
