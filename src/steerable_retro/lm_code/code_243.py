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
    This function detects a strategy where a flexible ethylene linker is converted
    to a rigid alkyne, creating conformational constraint in the molecule.
    """
    has_alkyne_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_alkyne_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for alkyne formation
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    product_mol = Chem.MolFromSmiles(product)

                    if reactant_mol and product_mol:
                        # Check if product has alkyne but reactant doesn't
                        alkyne_pattern = Chem.MolFromSmarts("[C]#[C]")
                        ethyl_pattern = Chem.MolFromSmarts("[CH2][CH2]")

                        if (
                            product_mol.HasSubstructMatch(alkyne_pattern)
                            and not reactant_mol.HasSubstructMatch(alkyne_pattern)
                            and reactant_mol.HasSubstructMatch(ethyl_pattern)
                        ):
                            print(f"Found alkyne formation at depth {depth}")
                            has_alkyne_formation = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return has_alkyne_formation
