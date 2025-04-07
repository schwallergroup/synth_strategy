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
    This function detects if the synthesis uses a biaryl formation strategy
    via coupling reactions (like Suzuki coupling).
    """
    biaryl_formation_detected = False

    def dfs_traverse(node):
        nonlocal biaryl_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for potential coupling reaction patterns
                # Look for boronic acid/ester and halide in reactants
                has_boronic = False
                has_halide = False

                boronic_pattern = Chem.MolFromSmarts("[c][B]([O])[O]")
                halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(boronic_pattern):
                            has_boronic = True
                        if reactant_mol.HasSubstructMatch(halide_pattern):
                            has_halide = True

                # Check if product has a new biaryl bond
                if has_boronic and has_halide:
                    print(
                        "Potential coupling reaction detected with boronic acid and halide"
                    )
                    biaryl_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    if biaryl_formation_detected:
        print("Biaryl formation via coupling strategy detected")
    else:
        print("No biaryl formation via coupling strategy detected")

    return biaryl_formation_detected
