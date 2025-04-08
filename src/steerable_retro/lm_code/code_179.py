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
    This function detects if the synthesis route involves a late-stage amide coupling
    (in the final or penultimate step).
    """
    amide_formation_at_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_at_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for amide formation
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Look for carboxylic acid in reactants
            acid_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#8;H1]")
            has_acid = any(
                mol.HasSubstructMatch(acid_pattern) for mol in reactant_mols if mol is not None
            )

            # Look for amine in reactants
            amine_pattern = Chem.MolFromSmarts("[#7;H2]-[#6]")
            has_amine = any(
                mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols if mol is not None
            )

            # Look for amide in product
            amide_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#7]-[#6]")
            has_amide_product = product_mol is not None and product_mol.HasSubstructMatch(
                amide_pattern
            )

            # If we have acid + amine → amide, it's an amide coupling
            if has_acid and has_amine and has_amide_product:
                amide_formation_at_depth = depth
                print(f"Detected amide coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if amide formation occurred at depth 0 or 1 (late stage)
    is_late_stage = amide_formation_at_depth is not None and amide_formation_at_depth <= 1
    print(f"Late stage amide coupling: {is_late_stage}")
    return is_late_stage
