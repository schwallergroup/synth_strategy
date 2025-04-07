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
    This function detects if the synthetic route contains a Suzuki coupling reaction
    for biaryl formation.
    """
    suzuki_detected = False

    def dfs_traverse(node):
        nonlocal suzuki_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl halide and boronic acid in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("[c]-[#53,#35,#17]")
            boronic_acid_pattern = Chem.MolFromSmarts("[#5]([#8])[#8]")

            # Check for biaryl formation in product
            biaryl_pattern = Chem.MolFromSmarts("c:c-c:c")

            has_aryl_halide = False
            has_boronic_acid = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True
                    if mol.HasSubstructMatch(boronic_acid_pattern):
                        has_boronic_acid = True

            product_mol = Chem.MolFromSmiles(product)
            has_biaryl = product_mol and product_mol.HasSubstructMatch(biaryl_pattern)

            if has_aryl_halide and has_boronic_acid and has_biaryl:
                print("Suzuki coupling detected")
                suzuki_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return suzuki_detected
