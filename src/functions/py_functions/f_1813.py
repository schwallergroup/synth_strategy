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
    Detects borylation of aryl halide as preparation for cross-coupling.
    Looks for conversion of aryl halide to boronate ester.
    """
    borylation_detected = False

    def dfs_traverse(node):
        nonlocal borylation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl halide in reactants
                aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")
                # Check for diboron reagent in reactants
                diboron_pattern = Chem.MolFromSmarts("[B][O][C].[B][O][C]")
                # Check for boronate in product
                boronate_pattern = Chem.MolFromSmarts("[c][B][O][C]")

                has_aryl_halide = False
                has_diboron = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True
                        if mol.HasSubstructMatch(diboron_pattern):
                            has_diboron = True

                prod_mol = Chem.MolFromSmiles(product)
                if (
                    prod_mol
                    and has_aryl_halide
                    and prod_mol.HasSubstructMatch(boronate_pattern)
                ):
                    borylation_detected = True
                    print(f"Detected borylation of aryl halide: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return borylation_detected
