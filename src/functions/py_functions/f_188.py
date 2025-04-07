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
    This function detects a late-stage amide formation involving a terminal alkyne.
    """
    has_late_stage_amide = False

    def dfs_traverse(node):
        nonlocal has_late_stage_amide

        if (
            node["type"] == "reaction" and node.get("depth", 0) == 0
        ):  # Depth 0 is the last step
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an amide formation
                amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                    # Check if one of the reactants has a terminal alkyne
                    terminal_alkyne_pattern = Chem.MolFromSmarts("[C]#[C][#6]")

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            terminal_alkyne_pattern
                        ):
                            has_late_stage_amide = True
                            print(
                                "Detected late-stage amide formation with terminal alkyne"
                            )
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_late_stage_amide
