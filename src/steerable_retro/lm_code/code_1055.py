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
    Detects if the synthesis route uses a Suzuki coupling strategy for biaryl formation.
    """
    found_suzuki = False

    def dfs_traverse(node, depth=0):
        nonlocal found_suzuki

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for patterns indicative of Suzuki coupling
                has_boronic = False
                has_halide = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Check for boronic acid/ester
                        if "B(O)" in reactant or "B(OC)" in reactant:
                            has_boronic = True

                        # Check for aryl halide
                        if mol.HasSubstructMatch(Chem.MolFromSmarts("c[Cl,Br,I]")):
                            has_halide = True

                # Check if product has a biaryl bond that wasn't in the reactants
                if has_boronic and has_halide:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts("c:c")):
                        found_suzuki = True
                        print(f"Found Suzuki coupling at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_suzuki
