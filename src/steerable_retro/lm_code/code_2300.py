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
    Detects if the route employs a Suzuki coupling strategy for C-C bond formation
    between an aryl halide and a boronic acid.
    """
    found_suzuki = False

    def dfs_traverse(node, depth=0):
        nonlocal found_suzuki

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl halide and boronic acid in reactants
                aryl_halide_pattern = Chem.MolFromSmarts("c[Br,I]")
                boronic_acid_pattern = Chem.MolFromSmarts("cB(O)O")

                has_aryl_halide = False
                has_boronic_acid = False

                for r in reactants:
                    try:
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            if mol.HasSubstructMatch(aryl_halide_pattern):
                                has_aryl_halide = True
                            if mol.HasSubstructMatch(boronic_acid_pattern):
                                has_boronic_acid = True
                    except:
                        continue

                if has_aryl_halide and has_boronic_acid:
                    found_suzuki = True
                    print(f"Found Suzuki coupling at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_suzuki
