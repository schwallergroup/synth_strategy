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
    Detects if the synthesis route maintains a tert-butyl ester protecting group
    throughout most of the synthesis.
    """
    tert_butyl_ester_count = 0
    reaction_count = 0

    def dfs_traverse(node):
        nonlocal tert_butyl_ester_count, reaction_count

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            reaction_count += 1
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Check for tert-butyl ester in product
            tert_butyl_ester_pattern = Chem.MolFromSmarts("C(=O)OC(C)(C)C")

            try:
                mol = Chem.MolFromSmiles(product)
                if mol and mol.HasSubstructMatch(tert_butyl_ester_pattern):
                    tert_butyl_ester_count += 1
                    print(f"Found tert-butyl ester in reaction product")
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Check if tert-butyl ester is present in at least 70% of reactions
    if reaction_count > 0 and tert_butyl_ester_count / reaction_count >= 0.7:
        print(
            f"tert-butyl ester protection strategy detected ({tert_butyl_ester_count}/{reaction_count} reactions)"
        )
        return True
    return False
