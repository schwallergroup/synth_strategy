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
    Detects if the synthetic route involves sequential functionalization:
    aryl halide → alkyne coupling → N-alkylation
    """
    # Track steps found
    found_aryl_halide = False
    found_alkyne_coupling = False
    found_n_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_aryl_halide, found_alkyne_coupling, found_n_alkylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl halide in early steps (higher depth)
                if depth >= 2:
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            aryl_halide_pattern = Chem.MolFromSmarts("[c]-[#53]")
                            if reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                                found_aryl_halide = True
                                print(f"Found aryl halide at depth {depth}")

                # Check for alkyne coupling in middle steps
                if 1 <= depth <= 2:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        alkyne_pattern = Chem.MolFromSmarts("[c]-[C]#[C]")
                        if product_mol.HasSubstructMatch(alkyne_pattern):
                            found_alkyne_coupling = True
                            print(f"Found alkyne coupling at depth {depth}")

                # Check for N-alkylation in late steps (lower depth)
                if depth <= 1:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        n_alkylation_pattern = Chem.MolFromSmarts("[n]([C])[c]")
                        piperidine_pattern = Chem.MolFromSmarts("[N]1[C][C][C][C][C]1")
                        if product_mol.HasSubstructMatch(
                            n_alkylation_pattern
                        ) and product_mol.HasSubstructMatch(piperidine_pattern):
                            found_n_alkylation = True
                            print(f"Found N-alkylation at depth {depth}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True only if all three steps are found in the correct sequence
    return found_aryl_halide and found_alkyne_coupling and found_n_alkylation
