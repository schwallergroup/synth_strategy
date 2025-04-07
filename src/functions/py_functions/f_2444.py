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
    Detects if the synthesis route contains a biaryl formation via Suzuki coupling
    (reaction between aryl halide and boronic acid).
    """
    found = False

    def dfs_traverse(node, depth=0):
        nonlocal found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl halide and boronic acid in reactants
                has_aryl_halide = False
                has_boronic_acid = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Check for aryl halide
                        if mol.HasSubstructMatch(Chem.MolFromSmarts("[c][Br,I]")):
                            has_aryl_halide = True
                        # Check for boronic acid
                        if mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[c][B]([O][H])[O][H]")
                        ) or mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[c][B]([O])[O]")
                        ):
                            has_boronic_acid = True

                # Check if product has a biaryl system not present in reactants
                if has_aryl_halide and has_boronic_acid:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[c]!@[c]")
                    ):
                        found = True
                        print(
                            f"Found biaryl formation via Suzuki coupling at depth {depth}"
                        )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found
