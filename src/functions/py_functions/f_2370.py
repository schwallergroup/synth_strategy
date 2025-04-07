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
    This function detects a convergent synthesis strategy using Suzuki coupling.
    It looks for reactions where a boronic acid/ester and aryl halide are combined.
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

            # Check if we have multiple reactants (convergent)
            if len(reactants) > 1:
                # Look for boronic acid/ester pattern in any reactant
                boronic_pattern = Chem.MolFromSmarts("[#6]-[#5;X3,X4]")
                # Look for aryl halide pattern in any reactant
                halide_pattern = Chem.MolFromSmarts("[#6]:[#6]-[#35,#53,#17]")

                has_boronic = False
                has_halide = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(boronic_pattern):
                                has_boronic = True
                            if mol.HasSubstructMatch(halide_pattern):
                                has_halide = True
                    except:
                        continue

                # If both patterns are found in different reactants, it's likely a Suzuki coupling
                if has_boronic and has_halide:
                    print("Suzuki coupling detected in reaction with RSMI:", rsmi)
                    suzuki_detected = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return suzuki_detected
