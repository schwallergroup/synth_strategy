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
    Detects if the synthesis involves a Suzuki coupling reaction forming
    a C-C bond between two aromatic systems.
    """
    suzuki_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Look for boronic acid/ester pattern in reactants
                boronic_pattern = Chem.MolFromSmarts("[#6]-[#5](-[#8])-[#8]")

                # Look for aryl halide pattern in reactants
                aryl_halide_pattern = Chem.MolFromSmarts("c-[#9,#17,#35,#53]")

                has_boronic = False
                has_aryl_halide = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(boronic_pattern):
                            has_boronic = True
                            print(
                                f"Found boronic acid/ester in reactant at depth {depth}"
                            )
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True
                            print(f"Found aryl halide in reactant at depth {depth}")

                # If both patterns are found in reactants, likely a Suzuki coupling
                if has_boronic and has_aryl_halide:
                    suzuki_detected = True
                    print(f"Suzuki coupling detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return suzuki_detected
