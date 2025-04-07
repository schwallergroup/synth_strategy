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
    (aryl halide + boronate ester → biaryl).
    """
    suzuki_coupling_detected = False

    def dfs_traverse(node):
        nonlocal suzuki_coupling_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Suzuki coupling pattern
            boronate_patt = Chem.MolFromSmarts("[#6]-[B]([O][C])[O][C]")
            halide_patt = Chem.MolFromSmarts("[#6]-[Br,I,Cl]")

            has_boronate = False
            has_halide = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if not reactant_mol:
                    continue

                if reactant_mol.HasSubstructMatch(boronate_patt):
                    has_boronate = True
                if reactant_mol.HasSubstructMatch(halide_patt):
                    has_halide = True

            if has_boronate and has_halide:
                # Check if product has a new biaryl bond
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    print("Suzuki coupling detected")
                    suzuki_coupling_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_coupling_detected
