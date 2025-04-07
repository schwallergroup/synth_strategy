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
    Detects if the synthetic route involves multiple Sonogashira couplings
    (C-C bond formation with alkynes) in sequence.
    """
    sonogashira_count = 0

    def dfs_traverse(node):
        nonlocal sonogashira_count

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Sonogashira coupling pattern
            # Look for formation of aryl-alkyne bond
            product_mol = Chem.MolFromSmiles(product)

            # Check if any reactant has terminal alkyne or TMS-alkyne
            has_alkyne = False
            for reactant in reactants:
                r_mol = Chem.MolFromSmiles(reactant)
                if r_mol:
                    # Terminal alkyne pattern
                    if r_mol.HasSubstructMatch(Chem.MolFromSmarts("[C]#[CH]")):
                        has_alkyne = True
                    # TMS-alkyne pattern
                    elif r_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[C]#[C]-[Si]([C])([C])[C]")
                    ):
                        has_alkyne = True

            # Check if any reactant has aryl-halide
            has_aryl_halide = False
            for reactant in reactants:
                r_mol = Chem.MolFromSmiles(reactant)
                if r_mol:
                    if r_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[Br,I,Cl]")):
                        has_aryl_halide = True

            # Check if product has aryl-alkyne bond
            if product_mol and has_alkyne and has_aryl_halide:
                if product_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[C]#[C]")):
                    sonogashira_count += 1
                    print(
                        f"Detected Sonogashira coupling at depth {node.get('depth', 'unknown')}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return sonogashira_count >= 2
