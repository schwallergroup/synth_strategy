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
    This function detects functional groups (cyano and fluorine) that are maintained
    throughout the synthesis from early to final stages.
    """
    final_product_has_cyano = False
    final_product_has_fluorine = False
    early_stage_has_cyano = False
    early_stage_has_fluorine = False

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_cyano, final_product_has_fluorine
        nonlocal early_stage_has_cyano, early_stage_has_fluorine

        if node["type"] == "mol" and node["smiles"]:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for cyano group
                cyano_pattern = Chem.MolFromSmarts("C#N")
                has_cyano = mol.HasSubstructMatch(cyano_pattern)

                # Check for fluorine
                fluorine_pattern = Chem.MolFromSmarts("[F]")
                has_fluorine = mol.HasSubstructMatch(fluorine_pattern)

                # Final product (depth 0)
                if depth == 0:
                    if has_cyano:
                        final_product_has_cyano = True
                        print(f"Final product has cyano group: {node['smiles']}")
                    if has_fluorine:
                        final_product_has_fluorine = True
                        print(f"Final product has fluorine: {node['smiles']}")

                # Early stage (depth >= 3)
                if depth >= 3:
                    if has_cyano:
                        early_stage_has_cyano = True
                        print(f"Early stage intermediate has cyano group: {node['smiles']}")
                    if has_fluorine:
                        early_stage_has_fluorine = True
                        print(f"Early stage intermediate has fluorine: {node['smiles']}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Return True if both functional groups are maintained throughout
    return (
        final_product_has_cyano
        and early_stage_has_cyano
        and final_product_has_fluorine
        and early_stage_has_fluorine
    )
