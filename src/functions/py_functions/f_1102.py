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
    This function detects if the synthesis includes a late-stage N-cyanation.
    """
    late_stage_cyanation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_cyanation_found

        if node["type"] == "reaction" and depth <= 1:  # Late in synthesis (lower depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a cyanation reaction
                cyano_pattern = Chem.MolFromSmarts("[N]C#N")
                bromo_cyano_pattern = Chem.MolFromSmarts("BrC#N")

                # Check if cyano group in product and bromo-cyano in reactants
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(cyano_pattern):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            bromo_cyano_pattern
                        ):
                            late_stage_cyanation_found = True
                            print(f"Late-stage cyanation found at depth {depth}")
                            break

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return late_stage_cyanation_found
