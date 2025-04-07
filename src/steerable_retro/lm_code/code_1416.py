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
    This function detects if the synthetic route involves late-stage N-alkylation of an indole.
    """
    indole_pattern = Chem.MolFromSmarts("c1cccc2[nH]ccc12")
    n_alkylated_indole_pattern = Chem.MolFromSmarts("c1cccc2n(C)ccc12")

    late_stage_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_alkylation

        if node["type"] == "reaction" and depth <= 1:  # Focus on late-stage reactions (low depth)
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an N-alkylation of indole
            has_indole_reactant = any(
                Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).HasSubstructMatch(indole_pattern)
                for r in reactants
            )
            has_alkylated_product = Chem.MolFromSmiles(product) and Chem.MolFromSmiles(
                product
            ).HasSubstructMatch(n_alkylated_indole_pattern)

            if has_indole_reactant and has_alkylated_product:
                late_stage_alkylation = True
                print(f"Late-stage indole N-alkylation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return late_stage_alkylation
