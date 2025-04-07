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
    This function detects if the synthetic route uses reductive amination
    as a late-stage coupling strategy (at depth 0-1).
    """
    late_stage_reductive_amination_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_reductive_amination_found

        if node["type"] == "reaction" and depth <= 1:  # Only check at depths 0-1
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for reductive amination pattern
            # Look for aldehyde in reactants
            aldehyde_pattern = Chem.MolFromSmarts("[#6]=O")
            # Look for amine in reactants
            amine_pattern = Chem.MolFromSmarts("[#7;H2]")
            # Look for secondary amine in product
            sec_amine_pattern = Chem.MolFromSmarts("[#6]-[#7;H1]-[#6]")

            if (
                product
                and any(r for r in reactants if r and r.HasSubstructMatch(aldehyde_pattern))
                and any(r for r in reactants if r and r.HasSubstructMatch(amine_pattern))
                and product.HasSubstructMatch(sec_amine_pattern)
            ):
                late_stage_reductive_amination_found = True
                print(f"Late-stage reductive amination detected at depth {depth}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Late-stage reductive amination strategy detected: {late_stage_reductive_amination_found}"
    )
    return late_stage_reductive_amination_found
