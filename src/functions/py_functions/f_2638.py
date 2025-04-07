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
    Detects if the synthesis uses Sonogashira coupling with TMS-protected alkyne.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # Check for Sonogashira coupling pattern
                if len(reactants) >= 2:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                    product_mol = Chem.MolFromSmiles(product)

                    # Patterns to look for
                    aryl_halide_pattern = Chem.MolFromSmarts("[c][I,Br,Cl]")
                    tms_alkyne_pattern = Chem.MolFromSmarts("[C][Si]([C])([C])[C]#[C]")

                    has_aryl_halide = False
                    has_tms_alkyne = False

                    for r_mol in reactant_mols:
                        if r_mol:
                            if r_mol.HasSubstructMatch(aryl_halide_pattern):
                                has_aryl_halide = True
                            if r_mol.HasSubstructMatch(tms_alkyne_pattern):
                                has_tms_alkyne = True

                    if has_aryl_halide and has_tms_alkyne:
                        # Check if product has aryl-alkyne connection
                        aryl_alkyne_pattern = Chem.MolFromSmarts("[c][C]#[C]")
                        if product_mol and product_mol.HasSubstructMatch(
                            aryl_alkyne_pattern
                        ):
                            print("Detected Sonogashira coupling with TMS-alkyne")
                            result = True
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return result
