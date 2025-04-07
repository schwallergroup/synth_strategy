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
    Detects if the synthetic route employs a Sonogashira coupling strategy
    (aryl halide + terminal alkyne → aryl-alkyne)
    """
    sonogashira_detected = False

    def dfs_traverse(node):
        nonlocal sonogashira_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl iodide in reactants
            aryl_iodide_pattern = Chem.MolFromSmarts("[c]-[I]")
            # Check for terminal alkyne in reactants
            terminal_alkyne_pattern = Chem.MolFromSmarts("[C]#[CH]")
            # Check for aryl-alkyne in product
            aryl_alkyne_pattern = Chem.MolFromSmarts("[c]-[C]#[C]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if (
                product_mol
                and any(
                    mol and mol.HasSubstructMatch(aryl_iodide_pattern)
                    for mol in reactant_mols
                )
                and any(
                    mol and mol.HasSubstructMatch(terminal_alkyne_pattern)
                    for mol in reactant_mols
                )
                and product_mol.HasSubstructMatch(aryl_alkyne_pattern)
            ):
                print("Sonogashira coupling detected")
                sonogashira_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return sonogashira_detected
