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
    Detects if the synthetic route employs a Suzuki coupling to form a biaryl system.
    """
    has_suzuki_coupling = False

    def dfs_traverse(node):
        nonlocal has_suzuki_coupling

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Split into reactants and product
                parts = rsmi.split(">")
                if len(parts) >= 3:
                    reactants = parts[0].split(".")
                    product = parts[-1]

                    # Check for boronic acid in reactants
                    boronic_acid_pattern = Chem.MolFromSmarts("[#6]B(O)O")

                    # Check for aryl halide in reactants
                    aryl_halide_pattern = Chem.MolFromSmarts("[#6]!@[#53,#35,#17]")

                    has_boronic_acid = False
                    has_aryl_halide = False

                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                if mol.HasSubstructMatch(boronic_acid_pattern):
                                    has_boronic_acid = True
                                if mol.HasSubstructMatch(aryl_halide_pattern):
                                    has_aryl_halide = True
                        except:
                            continue

                    # Check if product has a biaryl system
                    if has_boronic_acid and has_aryl_halide:
                        try:
                            product_mol = Chem.MolFromSmiles(product)
                            # Simple check for biaryl - two aromatic rings connected
                            biaryl_pattern = Chem.MolFromSmarts("c!@c")
                            if product_mol and product_mol.HasSubstructMatch(
                                biaryl_pattern
                            ):
                                print("Found Suzuki coupling forming biaryl system")
                                has_suzuki_coupling = True
                        except:
                            pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_suzuki_coupling
