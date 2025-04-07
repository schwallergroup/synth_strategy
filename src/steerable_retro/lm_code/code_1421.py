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
    This function detects a strategy involving late-stage functional group modification.
    """
    has_late_stage_modification = False

    def dfs_traverse(node):
        nonlocal has_late_stage_modification

        if node["type"] == "reaction":
            depth = node.get("depth", 0)

            # Consider depth 0 or 1 as late-stage
            if depth <= 1:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                # Check for functional group patterns
                functional_group_patterns = [
                    Chem.MolFromSmarts("O[CH2]c1ccccc1"),  # Benzyl ether
                    Chem.MolFromSmarts("[OH]c"),  # Phenol
                    Chem.MolFromSmarts("S(=O)(=O)[NH]"),  # Sulfonamide
                ]

                # Check if any functional group pattern is in reactants
                has_functional_group_reactant = any(
                    mol is not None
                    and any(mol.HasSubstructMatch(pattern) for pattern in functional_group_patterns)
                    for mol in reactant_mols
                )

                # Check if product has different functional groups than reactants
                if has_functional_group_reactant and product_mol is not None:
                    has_late_stage_modification = True
                    print(f"Found late-stage functional group modification at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_late_stage_modification
