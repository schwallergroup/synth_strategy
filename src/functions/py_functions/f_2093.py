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
    This function detects if the synthesis incorporates a thiazolidinedione-containing fragment,
    particularly in late-stage synthesis.
    """
    thiazolidinedione_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal thiazolidinedione_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for thiazolidinedione in reactants or product
            thiazolidinedione_pattern = Chem.MolFromSmarts("C1SC(=O)NC1=O")

            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol and product_mol.HasSubstructMatch(thiazolidinedione_pattern):
                # Check if it's newly incorporated
                has_thiazolidinedione_in_reactants = False
                for r_smiles in reactants_smiles:
                    r_mol = Chem.MolFromSmiles(r_smiles)
                    if r_mol and r_mol.HasSubstructMatch(thiazolidinedione_pattern):
                        has_thiazolidinedione_in_reactants = True
                        break

                if has_thiazolidinedione_in_reactants:
                    thiazolidinedione_depth = depth
                    print(f"Thiazolidinedione incorporation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if thiazolidinedione is incorporated in late-stage (low depth)
    is_late_stage_incorporation = (
        thiazolidinedione_depth >= 0 and thiazolidinedione_depth <= 1
    )

    print(
        f"Late-stage thiazolidinedione incorporation detected: {is_late_stage_incorporation}"
    )
    return is_late_stage_incorporation
