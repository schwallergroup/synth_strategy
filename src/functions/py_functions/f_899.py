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
    This function detects a sequence of Suzuki coupling reactions in a synthetic route.
    It looks for reactions where a boronic acid/ester and aryl halide form a biaryl product.
    """
    suzuki_count = 0

    def dfs_traverse(node):
        nonlocal suzuki_count

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for boronic acid/ester pattern in reactants
            boronic_pattern = Chem.MolFromSmarts("[c][B]([O])[O]")
            boronic_ester_pattern = Chem.MolFromSmarts(
                "[c][B]1[O][C]([C])([C])[C]([C])([C])[O]1"
            )

            # Check for aryl halide pattern in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")

            # Check for biaryl formation in product
            biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")

            has_boronic = False
            has_aryl_halide = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(
                            boronic_pattern
                        ) or mol.HasSubstructMatch(boronic_ester_pattern):
                            has_boronic = True
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product)
                has_biaryl = product_mol and product_mol.HasSubstructMatch(
                    biaryl_pattern
                )
            except:
                has_biaryl = False

            if has_boronic and has_aryl_halide and has_biaryl:
                print(f"Found Suzuki coupling: {rsmi}")
                suzuki_count += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Total Suzuki couplings found: {suzuki_count}")
    return suzuki_count >= 2  # Return True if at least 2 Suzuki couplings are found
