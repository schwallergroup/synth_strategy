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
    Detects if the synthesis route involves a sequence of transformations from
    nitro group to amine to thioamide to heterocycle.
    """
    # Track functional groups at each depth
    functional_groups_by_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for functional groups
                nitro_pattern = Chem.MolFromSmarts("[#8-][#7+](=[#8])[#6]")
                amine_pattern = Chem.MolFromSmarts("[NH2][#6]")
                thioamide_pattern = Chem.MolFromSmarts("[#7][#6](=[#16])[#7]")
                thiazole_pattern = Chem.MolFromSmarts("[#6]1[#16][#6][#6][#7]1")

                if depth not in functional_groups_by_depth:
                    functional_groups_by_depth[depth] = set()

                if mol.HasSubstructMatch(nitro_pattern):
                    functional_groups_by_depth[depth].add("nitro")
                if mol.HasSubstructMatch(amine_pattern):
                    functional_groups_by_depth[depth].add("amine")
                if mol.HasSubstructMatch(thioamide_pattern):
                    functional_groups_by_depth[depth].add("thioamide")
                if mol.HasSubstructMatch(thiazole_pattern):
                    functional_groups_by_depth[depth].add("thiazole")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check for the transformation sequence
    depths = sorted(functional_groups_by_depth.keys())

    # Check if we have the sequence nitro -> amine -> thioamide -> thiazole
    # Note: In retrosynthetic direction (increasing depth), we'd see thiazole -> thioamide -> amine -> nitro
    has_nitro = False
    has_amine = False
    has_thioamide = False
    has_thiazole = False

    for depth in depths:
        if "nitro" in functional_groups_by_depth[depth]:
            has_nitro = True
        if "amine" in functional_groups_by_depth[depth]:
            has_amine = True
        if "thioamide" in functional_groups_by_depth[depth]:
            has_thioamide = True
        if "thiazole" in functional_groups_by_depth[depth]:
            has_thiazole = True

    sequence_found = has_nitro and has_amine and has_thioamide and has_thiazole

    if sequence_found:
        print("Detected nitro->amine->thioamide->heterocycle transformation sequence")

    return sequence_found
