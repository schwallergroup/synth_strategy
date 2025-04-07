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
    This function detects Sonogashira coupling (aryl halide + terminal alkyne).
    """
    sonogashira_found = False

    def dfs_traverse(node):
        nonlocal sonogashira_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Sonogashira coupling
            terminal_alkyne_pattern = Chem.MolFromSmarts("[C]#[CH]")
            aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Cl,Br,I]")
            internal_alkyne_pattern = Chem.MolFromSmarts("[c]-[C]#[C]-[c]")

            try:
                # Check if product has internal alkyne
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(
                    internal_alkyne_pattern
                ):
                    # Check if reactants have terminal alkyne and aryl halide
                    terminal_alkyne_found = False
                    aryl_halide_found = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(terminal_alkyne_pattern):
                                terminal_alkyne_found = True
                            if reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                                aryl_halide_found = True

                    if terminal_alkyne_found and aryl_halide_found:
                        sonogashira_found = True
                        print("Sonogashira coupling detected")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return sonogashira_found
