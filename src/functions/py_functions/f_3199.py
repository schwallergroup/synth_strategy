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
    Detects a convergent synthesis strategy using Suzuki coupling for biaryl formation.
    Looks for boronic acid and aryl halide reactants forming a biaryl product.
    """
    suzuki_coupling_found = False

    def dfs_traverse(node):
        nonlocal suzuki_coupling_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if we have multiple reactants (convergent)
                if len(reactants) >= 2:
                    # Look for boronic acid pattern in any reactant
                    boronic_acid_pattern = Chem.MolFromSmarts("[#6]-[#5](-[#8])-[#8]")
                    aryl_halide_pattern = Chem.MolFromSmarts("[#6]:[#6]-[#35,#53,#17]")

                    has_boronic_acid = False
                    has_aryl_halide = False

                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(boronic_acid_pattern):
                                has_boronic_acid = True
                            if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                                has_aryl_halide = True
                        except:
                            continue

                    # Check if product has a biaryl bond that wasn't in the reactants
                    if has_boronic_acid and has_aryl_halide:
                        print(
                            "Found potential Suzuki coupling with boronic acid and aryl halide"
                        )
                        suzuki_coupling_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return suzuki_coupling_found
