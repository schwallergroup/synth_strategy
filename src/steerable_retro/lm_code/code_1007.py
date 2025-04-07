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
    (aryl halide + boronic acid → biaryl).
    """
    has_suzuki_coupling = False

    def dfs_traverse(node):
        nonlocal has_suzuki_coupling

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # Check for Suzuki coupling pattern
                product_mol = Chem.MolFromSmiles(product)

                if product_mol:
                    # Check if product contains biaryl bond
                    biaryl_pattern = Chem.MolFromSmarts("c:c-c:c")

                    if product_mol.HasSubstructMatch(biaryl_pattern):
                        # Check if reactants contain aryl halide and boronic acid
                        has_aryl_halide = False
                        has_boronic_acid = False

                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                aryl_halide_pattern = Chem.MolFromSmarts("c[I,Br,Cl]")
                                boronic_acid_pattern = Chem.MolFromSmarts("cB(O)O")

                                if reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                                    has_aryl_halide = True
                                if reactant_mol.HasSubstructMatch(boronic_acid_pattern):
                                    has_boronic_acid = True

                        if has_aryl_halide and has_boronic_acid:
                            has_suzuki_coupling = True
                            print(
                                f"Detected Suzuki coupling at depth {node.get('depth', 'unknown')}"
                            )
            except Exception as e:
                print(f"Error in Suzuki coupling detection: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return has_suzuki_coupling
