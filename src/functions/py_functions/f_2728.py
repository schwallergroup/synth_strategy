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
    Detects if the route uses a late-stage cross-coupling strategy (Suzuki coupling)
    early in the synthesis (high depth in retrosynthetic tree).
    """
    found_suzuki = False

    def dfs_traverse(node, depth=0):
        nonlocal found_suzuki

        if (
            node["type"] == "reaction" and depth >= 4
        ):  # Looking for cross-coupling at depth >= 4
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for boronic acid/ester pattern in reactants
                boronic_pattern = Chem.MolFromSmarts("[#6]-[B]([O,c,C])[O,c,C]")

                # Check for aryl halide pattern
                aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Br,I,Cl]")

                has_boronic = False
                has_aryl_halide = False

                for reactant in reactants:
                    if reactant.strip():
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(boronic_pattern):
                                has_boronic = True
                                print(f"Found boronic acid/ester at depth {depth}")
                            if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                                has_aryl_halide = True
                                print(f"Found aryl halide at depth {depth}")
                        except:
                            continue

                # Check if product has a new biaryl bond
                if has_boronic and has_aryl_halide:
                    try:
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol:
                            # This is a simplification - a more robust check would compare
                            # reactants and products to confirm new C-C bond formation
                            found_suzuki = True
                            print(f"Detected Suzuki coupling at depth {depth}")
                    except:
                        pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_suzuki
