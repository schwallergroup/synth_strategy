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
    This function detects late-stage reductive amination to form secondary amines.
    """
    found_reductive_amination = False
    max_depth = 0
    reductive_amination_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal found_reductive_amination, max_depth, reductive_amination_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for reductive amination (aldehyde + primary amine → secondary amine)
            aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
            primary_amine_pattern = Chem.MolFromSmarts("[NH2]")
            secondary_amine_pattern = Chem.MolFromSmarts("[NH]")

            has_aldehyde = any(
                Chem.MolFromSmiles(r)
                and Chem.MolFromSmiles(r).HasSubstructMatch(aldehyde_pattern)
                for r in reactants
                if r
            )
            has_primary_amine = any(
                Chem.MolFromSmiles(r)
                and Chem.MolFromSmiles(r).HasSubstructMatch(primary_amine_pattern)
                for r in reactants
                if r
            )
            product_mol = Chem.MolFromSmiles(product) if product else None
            has_secondary_amine = product_mol and product_mol.HasSubstructMatch(
                secondary_amine_pattern
            )

            if has_aldehyde and has_primary_amine and has_secondary_amine:
                found_reductive_amination = True
                reductive_amination_depth = depth
                print(f"Found reductive amination at depth {depth}")

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if reductive amination occurs late in the synthesis (in the first half of depth)
    if (
        found_reductive_amination
        and max_depth > 0
        and reductive_amination_depth is not None
    ):
        late_threshold = max_depth / 2
        late_stage = reductive_amination_depth < late_threshold
        return late_stage

    return False
