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
    This function detects a synthetic strategy involving late-stage reductive amination
    for C-N bond formation.
    """
    has_reductive_amination = False
    is_late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal has_reductive_amination, is_late_stage

        if node["type"] == "reaction":
            # Check if this is a low-depth (late in synthesis) reaction
            if depth <= 1:  # Depth 0 or 1 is considered late in synthesis
                is_late_stage = True

                # Extract reactants and product
                rsmi = node["metadata"].get("rsmi", "")
                if rsmi:
                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]

                    reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
                    product = Chem.MolFromSmiles(product_part) if product_part else None

                    if product and all(r for r in reactants):
                        # Check for aldehyde in reactants
                        aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
                        has_aldehyde = any(r.HasSubstructMatch(aldehyde_pattern) for r in reactants)

                        # Check for amine in reactants
                        amine_pattern = Chem.MolFromSmarts("[#7;H1]")
                        has_amine = any(r.HasSubstructMatch(amine_pattern) for r in reactants)

                        # Check if product has new C-N bond
                        if has_aldehyde and has_amine:
                            print("Detected potential reductive amination at depth", depth)
                            has_reductive_amination = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    result = has_reductive_amination and is_late_stage
    print(f"Late-stage reductive amination strategy detected: {result}")
    return result
