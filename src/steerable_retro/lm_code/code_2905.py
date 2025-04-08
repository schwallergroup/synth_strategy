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
    Detects if the synthetic route involves a Sonogashira-type coupling between
    an aryl halide and a terminal alkyne.
    """
    # Track if we found the pattern
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl halide in reactants
                aryl_halide_present = False
                terminal_alkyne_present = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        # Check for aryl iodide
                        aryl_iodide_pattern = Chem.MolFromSmarts("[c]-[#53]")
                        if reactant_mol.HasSubstructMatch(aryl_iodide_pattern):
                            aryl_halide_present = True

                        # Check for terminal alkyne (propargyl alcohol)
                        terminal_alkyne_pattern = Chem.MolFromSmarts("[C]#[C][C][O]")
                        if reactant_mol.HasSubstructMatch(terminal_alkyne_pattern):
                            terminal_alkyne_present = True

                # Check if product has C-C≡C where aryl-I was
                if aryl_halide_present and terminal_alkyne_present:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        alkyne_connection_pattern = Chem.MolFromSmarts("[c]-[C]#[C]")
                        if product_mol.HasSubstructMatch(alkyne_connection_pattern):
                            print("Found Sonogashira-type coupling with terminal alkyne")
                            found_pattern = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_pattern
