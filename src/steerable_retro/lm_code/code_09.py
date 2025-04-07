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
    Detects if the synthetic route involves SNAr with alcohol nucleophile
    to form aryl ether bonds.
    """
    snar_with_alcohol_found = False

    def dfs_traverse(node):
        nonlocal snar_with_alcohol_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for SNAr with alcohol
                has_aryl_halide = False
                has_alcohol = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        # Check for electron-deficient aryl halide (like pyridine-Cl)
                        if reactant_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[n]1[c]([Cl,Br,F,I])[c][c][c][c]1")
                        ) or reactant_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[c]1([Cl,Br,F,I])[n][c][c][c][c]1")
                        ):
                            has_aryl_halide = True

                        # Check for alcohol
                        if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[O][H]")):
                            has_alcohol = True

                # Check if product has aryl ether
                product_mol = Chem.MolFromSmiles(product)
                if (
                    has_aryl_halide
                    and has_alcohol
                    and product_mol
                    and product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[n]1[c]([O][#6])[c][c][c][c]1")
                    )
                    or product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[c]1([O][#6])[n][c][c][c][c]1")
                    )
                ):
                    print("Found SNAr with alcohol")
                    snar_with_alcohol_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return snar_with_alcohol_found
