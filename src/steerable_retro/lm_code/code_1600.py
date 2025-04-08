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
    Detects a synthetic strategy involving early-stage carbamate formation
    between aromatic components.
    """
    # Track if we found carbamate formation and at what depth
    carbamate_formation_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal carbamate_formation_depth, max_depth

        # Track maximum depth to determine early vs late stage
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if not product_mol or not all(reactant_mols):
                print("Warning: Could not parse some molecules in reaction")
                return

            # Check for carbamate formation
            phenol_pattern = Chem.MolFromSmarts("c[OH]")
            chloroformate_pattern = Chem.MolFromSmarts("ClC(=O)[O,N]")
            carbamate_pattern = Chem.MolFromSmarts("NC(=O)Oc")  # N-C(=O)-O-aryl pattern

            phenol_present = any(mol.HasSubstructMatch(phenol_pattern) for mol in reactant_mols)
            chloroformate_present = any(
                mol.HasSubstructMatch(chloroformate_pattern) for mol in reactant_mols
            )

            if (phenol_present or chloroformate_present) and product_mol.HasSubstructMatch(
                carbamate_pattern
            ):
                carbamate_formation_depth = depth
                print(f"Found carbamate formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we found carbamate formation and if it's in the early stage (high depth)
    early_stage_threshold = max_depth * 0.6  # Consider early stage if in the top 60% of depth
    strategy_present = (
        carbamate_formation_depth is not None and carbamate_formation_depth >= early_stage_threshold
    )

    print(
        f"Strategy detection result: {strategy_present} (depth: {carbamate_formation_depth}, threshold: {early_stage_threshold})"
    )
    return strategy_present
