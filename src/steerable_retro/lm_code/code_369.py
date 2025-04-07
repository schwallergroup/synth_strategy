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
    This function detects if the synthesis route involves a Strecker reaction
    (formation of α-aminonitrile from aldehyde, HCN, and ammonia).
    """
    strecker_found = False

    def dfs_traverse(node, depth=0):
        nonlocal strecker_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aldehyde in reactants
                aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
                # Check for HCN in reactants
                hcn_pattern = Chem.MolFromSmarts("[C]#[N]")
                # Check for ammonia in reactants
                ammonia_pattern = Chem.MolFromSmarts("[NH3]")
                # Check for α-aminonitrile in product
                aminonitrile_pattern = Chem.MolFromSmarts("[C]([NH2])[C]#[N]")

                aldehyde_present = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(aldehyde_pattern)
                    for r in reactants
                    if r
                )
                hcn_present = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(hcn_pattern)
                    for r in reactants
                    if r
                )
                ammonia_present = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(ammonia_pattern)
                    for r in reactants
                    if r
                )
                aminonitrile_in_product = Chem.MolFromSmiles(
                    product
                ) is not None and Chem.MolFromSmiles(product).HasSubstructMatch(
                    aminonitrile_pattern
                )

                if (
                    aldehyde_present
                    and (hcn_present or ammonia_present)
                    and aminonitrile_in_product
                ):
                    print(f"Found Strecker reaction at depth {depth}")
                    strecker_found = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return strecker_found
