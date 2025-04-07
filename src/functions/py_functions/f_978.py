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
    Detects if the synthesis route uses multiple Sonogashira couplings
    (C-C bond formation between terminal alkyne and aryl/heteroaryl halide).
    """
    sonogashira_count = 0

    def dfs_traverse(node):
        nonlocal sonogashira_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for terminal alkyne in reactants
                terminal_alkyne_pattern = Chem.MolFromSmarts("[C]#[CH]")

                # Check for aryl/heteroaryl halide in reactants
                aryl_halide_pattern = Chem.MolFromSmarts("c[I,Br,Cl]")
                heteroaryl_halide_pattern = Chem.MolFromSmarts("n[I,Br,Cl]")

                # Check reactants for patterns
                terminal_alkyne_present = False
                aryl_halide_present = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(terminal_alkyne_pattern):
                            terminal_alkyne_present = True
                        if mol.HasSubstructMatch(
                            aryl_halide_pattern
                        ) or mol.HasSubstructMatch(heteroaryl_halide_pattern):
                            aryl_halide_present = True

                # Check if product has new C-C bond between former alkyne and aryl
                if terminal_alkyne_present and aryl_halide_present:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol:
                        # If both patterns were in reactants and product has internal alkyne,
                        # likely a Sonogashira coupling occurred
                        internal_alkyne_pattern = Chem.MolFromSmarts("c[C]#[C]")
                        if prod_mol.HasSubstructMatch(internal_alkyne_pattern):
                            print("Sonogashira coupling detected")
                            sonogashira_count += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return (
        sonogashira_count >= 2
    )  # Return True if at least 2 Sonogashira couplings detected
