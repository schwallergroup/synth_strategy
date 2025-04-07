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
    This function detects if the synthetic route employs Suzuki coupling reactions
    (boronic acid/ester + aryl halide -> C-C bond formation)
    """
    suzuki_coupling_found = False

    def dfs_traverse(node):
        nonlocal suzuki_coupling_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if one reactant contains boronic acid/ester and another contains aryl halide
            has_boronic = False
            has_aryl_halide = False

            for reactant in reactants:
                if reactant:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Check for boronic acid/ester pattern
                        boronic_pattern = Chem.MolFromSmarts("[#6]-[#5](-[#8])-[#8]")
                        if mol.HasSubstructMatch(boronic_pattern):
                            has_boronic = True

                        # Check for aryl halide pattern
                        aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Br,I,Cl]")
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True

            # If both patterns are found in the reactants, it's likely a Suzuki coupling
            if has_boronic and has_aryl_halide:
                print("Suzuki coupling detected in reaction")
                suzuki_coupling_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return suzuki_coupling_found
