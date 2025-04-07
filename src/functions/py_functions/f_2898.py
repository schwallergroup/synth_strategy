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
    Detects if the synthetic route involves a Sonogashira coupling (alkyne + aryl halide).
    """
    sonogashira_found = False

    def dfs_traverse(node):
        nonlocal sonogashira_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Sonogashira coupling
                if len(reactants) >= 2:  # Need at least two reactants
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    product_mol = Chem.MolFromSmiles(product) if product else None

                    if (
                        all(m is not None for m in reactant_mols)
                        and product_mol is not None
                    ):
                        # Check for terminal alkyne
                        terminal_alkyne_pattern = Chem.MolFromSmarts("C#[CH]")
                        # Check for aryl halide
                        aryl_halide_pattern = Chem.MolFromSmarts("c-[#9,#17,#35,#53]")
                        # Check for coupled product
                        alkyne_aryl_pattern = Chem.MolFromSmarts("C#C-c")

                        has_terminal_alkyne = any(
                            len(m.GetSubstructMatches(terminal_alkyne_pattern)) > 0
                            for m in reactant_mols
                        )
                        has_aryl_halide = any(
                            len(m.GetSubstructMatches(aryl_halide_pattern)) > 0
                            for m in reactant_mols
                        )
                        has_coupled_product = product_mol.GetSubstructMatches(
                            alkyne_aryl_pattern
                        )

                        if (
                            has_terminal_alkyne
                            and has_aryl_halide
                            and has_coupled_product
                        ):
                            sonogashira_found = True
                            print("Found Sonogashira coupling reaction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return sonogashira_found
