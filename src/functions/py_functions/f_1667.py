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
    This function detects a synthetic strategy involving C-N bond formation
    between an aryl halide and an amine.
    """
    has_aryl_halide_coupling = False

    def dfs_traverse(node):
        nonlocal has_aryl_halide_coupling

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl halide and amine in reactants
            aryl_halide_present = False
            amine_present = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    aryl_halide_pattern = Chem.MolFromSmarts("[c][Cl,Br,I,F]")
                    amine_pattern = Chem.MolFromSmarts("[NH2][CH2]")

                    if mol.HasSubstructMatch(aryl_halide_pattern):
                        aryl_halide_present = True
                    if mol.HasSubstructMatch(amine_pattern):
                        amine_present = True

            # Check for C-N bond in product
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                c_n_bond_pattern = Chem.MolFromSmarts("[c][NH][CH2]")
                if product_mol.HasSubstructMatch(c_n_bond_pattern):
                    if aryl_halide_present and amine_present:
                        print("Found aryl halide coupling with amine")
                        has_aryl_halide_coupling = True

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_aryl_halide_coupling
