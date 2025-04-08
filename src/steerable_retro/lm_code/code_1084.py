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
    This function detects if the synthetic route uses a late-stage Sonogashira coupling
    (at depth 0 or 1) to connect major fragments.
    """
    has_late_stage_sonogashira = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_sonogashira

        if node["type"] == "reaction" and depth <= 1:
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
                            has_late_stage_sonogashira = True
                            print(
                                f"Detected late-stage Sonogashira coupling at depth {depth}: {rsmi}"
                            )
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_late_stage_sonogashira
