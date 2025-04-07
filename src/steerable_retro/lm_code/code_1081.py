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
    This function detects if the synthetic route uses Sonogashira coupling reactions
    to form C-C bonds between aryl halides and terminal alkynes.
    """
    sonogashira_count = 0

    def dfs_traverse(node):
        nonlocal sonogashira_count

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl halide and alkyne in reactants
                aryl_halide_pattern = Chem.MolFromSmarts("[c]-[#9,#17,#35,#53]")
                alkyne_pattern = Chem.MolFromSmarts("[C]#[C]")

                has_aryl_halide = False
                has_alkyne = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True
                        if mol and mol.HasSubstructMatch(alkyne_pattern):
                            has_alkyne = True
                    except:
                        continue

                # Check for aryl-alkyne in product
                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(alkyne_pattern):
                        if has_aryl_halide and has_alkyne:
                            sonogashira_count += 1
                            print(f"Detected Sonogashira coupling: {rsmi}")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return sonogashira_count >= 1
