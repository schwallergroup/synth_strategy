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
    This function detects if the synthesis route includes a Sonogashira-type coupling
    (C(sp²)-C(sp) bond formation between aryl halide and alkyne).
    """
    sonogashira_found = False

    def dfs_traverse(node):
        nonlocal sonogashira_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl halide and alkyne in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("[c]-[#9,#17,#35,#53]")
            alkyne_pattern = Chem.MolFromSmarts("[#6]#[#6]")

            aryl_halide_found = False
            alkyne_found = False

            try:
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                            aryl_halide_found = True
                        if reactant_mol.HasSubstructMatch(alkyne_pattern):
                            alkyne_found = True

                # Check if product has C(sp²)-C(sp) bond
                if aryl_halide_found and alkyne_found:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[c]-[#6]#[#6]")
                    ):
                        print("Sonogashira-type coupling detected")
                        sonogashira_found = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return sonogashira_found
