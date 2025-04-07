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
    This function detects if certain functional groups (like cyano) are preserved
    throughout the synthesis.
    """
    # Track functional groups at each depth
    functional_groups_by_depth = {}
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth

        if depth > max_depth:
            max_depth = depth

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])

            if mol:
                # Check for cyano group
                cyano_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")
                has_cyano = mol.HasSubstructMatch(cyano_pattern)

                # Check for halogens
                halogen_pattern = Chem.MolFromSmarts("[#6][Cl,Br,I]")
                has_halogen = mol.HasSubstructMatch(halogen_pattern)

                if depth not in functional_groups_by_depth:
                    functional_groups_by_depth[depth] = {
                        "cyano": has_cyano,
                        "halogen": has_halogen,
                    }
                else:
                    functional_groups_by_depth[depth]["cyano"] = (
                        functional_groups_by_depth[depth]["cyano"] or has_cyano
                    )
                    functional_groups_by_depth[depth]["halogen"] = (
                        functional_groups_by_depth[depth]["halogen"] or has_halogen
                    )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if cyano group is preserved throughout
    cyano_preserved = True
    halogen_preserved = True

    for depth in range(max_depth + 1):
        if depth in functional_groups_by_depth:
            if not functional_groups_by_depth[depth]["cyano"]:
                cyano_preserved = False
            if not functional_groups_by_depth[depth]["halogen"]:
                halogen_preserved = False

    if cyano_preserved:
        print("Cyano group preserved throughout synthesis")
    if halogen_preserved:
        print("Halogen groups preserved throughout synthesis")

    return cyano_preserved or halogen_preserved
