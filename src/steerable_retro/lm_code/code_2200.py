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
    This function detects a strategy involving the formation of a heterocycle
    (specifically a chroman ring) during the synthesis.
    """
    # Flag to track heterocycle formation
    has_heterocycle_formation = False

    def dfs_traverse(node):
        nonlocal has_heterocycle_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for heterocycle formation
                reactant_mol = Chem.MolFromSmiles(reactants[0]) if reactants else None
                product_mol = Chem.MolFromSmiles(product) if product else None

                if reactant_mol and product_mol:
                    # Check for chroman formation
                    chroman_pattern = Chem.MolFromSmarts("c1ccc2c(c1)CCCO2")

                    # Count rings in reactant and product
                    if not reactant_mol.HasSubstructMatch(
                        chroman_pattern
                    ) and product_mol.HasSubstructMatch(chroman_pattern):
                        print("Detected heterocycle (chroman) formation")
                        has_heterocycle_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_heterocycle_formation
