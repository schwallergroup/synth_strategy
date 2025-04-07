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
    Detects biaryl formation via Suzuki coupling (aryl halide + boronic acid/ester).
    """
    # Track if we found the pattern
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction":
            # Check if this is a reaction node
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Suzuki coupling patterns
                aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Br,I,Cl]")
                boronic_acid_pattern = Chem.MolFromSmarts("[c]-[B]([O])[O]")
                biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")

                # Check reactants for aryl halide and boronic acid
                has_aryl_halide = False
                has_boronic_acid = False

                for r in reactants:
                    try:
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            if mol.HasSubstructMatch(aryl_halide_pattern):
                                has_aryl_halide = True
                            if mol.HasSubstructMatch(boronic_acid_pattern):
                                has_boronic_acid = True
                    except:
                        continue

                # Check product for biaryl formation
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    has_biaryl = prod_mol and prod_mol.HasSubstructMatch(biaryl_pattern)
                except:
                    has_biaryl = False

                # If we have aryl halide, boronic acid, and biaryl in product, it's likely Suzuki coupling
                if has_aryl_halide and has_boronic_acid and has_biaryl:
                    found_pattern = True
                    print(
                        f"Found biaryl formation via Suzuki coupling at depth {depth}"
                    )

        # Traverse children
        if "children" in node:
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Biaryl formation via Suzuki coupling: {found_pattern}")
    return found_pattern
