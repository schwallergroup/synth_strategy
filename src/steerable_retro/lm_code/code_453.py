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
    Detects if the synthesis involves halogenation (I, Br) to prepare a coupling partner.
    """
    found_halogenation = False

    def dfs_traverse(node):
        nonlocal found_halogenation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]
                reactants = rsmi.split(">")[0].split(".")

                # Check if product contains aryl/heteroaryl halide
                product_mol = Chem.MolFromSmiles(product)
                aryl_halide_pattern = Chem.MolFromSmarts("c[Br,I]")

                if product_mol and product_mol.HasSubstructMatch(aryl_halide_pattern):
                    # Check if reactants don't have the halide
                    has_halide_in_reactants = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                            has_halide_in_reactants = True
                            break

                    if not has_halide_in_reactants:
                        print("Found halogenation to prepare coupling partner")
                        found_halogenation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_halogenation
