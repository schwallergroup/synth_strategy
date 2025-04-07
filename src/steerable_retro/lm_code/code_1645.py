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
    This function detects Suzuki coupling (aryl halide + boronic acid/ester).
    """
    has_suzuki_coupling = False

    def dfs_traverse(node):
        nonlocal has_suzuki_coupling

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl halide
                aryl_halide_pattern = Chem.MolFromSmarts("c[Br,I,Cl]")

                # Check for boronic acid/ester
                boronic_pattern = Chem.MolFromSmarts("[#6]B([#8])[#8]")
                boronic_ester_pattern = Chem.MolFromSmarts("[#6]B1OC(C)(C)OC1(C)C")

                has_aryl_halide = False
                has_boronic = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(aryl_halide_pattern):
                                has_aryl_halide = True
                            if mol.HasSubstructMatch(boronic_pattern) or mol.HasSubstructMatch(
                                boronic_ester_pattern
                            ):
                                has_boronic = True
                    except:
                        continue

                if has_aryl_halide and has_boronic:
                    print("Detected Suzuki coupling")
                    has_suzuki_coupling = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return has_suzuki_coupling
