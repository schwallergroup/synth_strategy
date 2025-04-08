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
    Detects if the synthesis route involves a Suzuki coupling reaction
    (boronic ester + aryl halide → C-C bond).
    """
    found_suzuki = False

    def dfs_traverse(node):
        nonlocal found_suzuki

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check for boronic ester pattern
            boronic_pattern = Chem.MolFromSmarts("[#5][#8][#6]")
            # Check for aryl halide pattern
            aryl_halide_pattern = Chem.MolFromSmarts("c[Br,I,Cl]")

            has_boronic = False
            has_aryl_halide = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if not reactant_mol:
                    continue

                if reactant_mol.HasSubstructMatch(boronic_pattern):
                    has_boronic = True

                if reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                    has_aryl_halide = True

            if has_boronic and has_aryl_halide:
                found_suzuki = True
                print(f"Found Suzuki coupling at depth: {node.get('depth', 'unknown')}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return found_suzuki
