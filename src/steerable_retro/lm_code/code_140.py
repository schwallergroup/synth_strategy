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
    Detects if the synthesis route includes a Suzuki coupling
    (formation of biaryl C-C bond from aryl halide and boronic acid)
    """
    found_suzuki = False

    def dfs_traverse(node):
        nonlocal found_suzuki

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Create RDKit mol objects
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                if product_mol and all(reactant_mols):
                    # SMARTS for aryl halide (Br, I, Cl)
                    aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")

                    # SMARTS for boronic acid or ester
                    boronic_pattern = Chem.MolFromSmarts("[c][B]([O])[O]")

                    # SMARTS for biaryl bond
                    biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")

                    # Check if reactants contain aryl halide and boronic acid
                    has_aryl_halide = any(
                        r.HasSubstructMatch(aryl_halide_pattern) for r in reactant_mols
                    )
                    has_boronic = any(r.HasSubstructMatch(boronic_pattern) for r in reactant_mols)

                    # Check if product has biaryl bond
                    has_biaryl = product_mol.HasSubstructMatch(biaryl_pattern)

                    if has_aryl_halide and has_boronic and has_biaryl:
                        print("Found Suzuki coupling")
                        found_suzuki = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return found_suzuki
