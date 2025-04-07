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
    This function detects if the synthetic route involves early TMS protection
    that is maintained through multiple steps of the synthesis.
    """
    tms_pattern = Chem.MolFromSmarts("[#6][Si]([#6])([#6])[#6]")
    tms_reactions = []
    tms_depth = -1
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal tms_reactions, tms_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if TMS is being introduced in this reaction
            if (
                not any(
                    Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).HasSubstructMatch(tms_pattern)
                    for r in reactants
                )
                and Chem.MolFromSmiles(product)
                and Chem.MolFromSmiles(product).HasSubstructMatch(tms_pattern)
            ):
                tms_reactions.append(depth)
                tms_depth = depth
                print(f"TMS introduction detected at depth {depth}")

            # Check if TMS is maintained in this reaction
            elif (
                any(
                    Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).HasSubstructMatch(tms_pattern)
                    for r in reactants
                )
                and Chem.MolFromSmiles(product)
                and Chem.MolFromSmiles(product).HasSubstructMatch(tms_pattern)
            ):
                tms_reactions.append(depth)
                print(f"TMS maintained at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if TMS was introduced early and maintained through multiple steps
    if tms_depth >= 0 and tms_depth >= max_depth - 1 and len(tms_reactions) >= 2:
        print(
            f"TMS protection strategy detected: introduced at depth {tms_depth} and maintained through {len(tms_reactions)} reactions"
        )
        return True
    return False
