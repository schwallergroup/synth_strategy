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
    This function detects if the synthesis incorporates a piperazine moiety in the late stage
    (second half) of the synthesis via SNAr or similar reaction.
    """
    piperazine_incorporation_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal piperazine_incorporation_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Check if this reaction incorporates a piperazine
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant contains piperazine
                piperazine_pattern = Chem.MolFromSmarts("[N]1CCN([*])CC1")
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(piperazine_pattern):
                            # Check if product also has piperazine (indicating incorporation)
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol and product_mol.HasSubstructMatch(
                                piperazine_pattern
                            ):
                                piperazine_incorporation_depth = depth
                                print(
                                    f"Piperazine incorporation detected at depth {depth}"
                                )
                    except:
                        continue

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Late stage is defined as occurring in the second half of the synthesis
    if piperazine_incorporation_depth is not None:
        is_late_stage = piperazine_incorporation_depth <= (max_depth / 2)
        print(
            f"Piperazine incorporation at depth {piperazine_incorporation_depth}, max depth {max_depth}"
        )
        print(f"Is late stage: {is_late_stage}")
        return is_late_stage

    return False
