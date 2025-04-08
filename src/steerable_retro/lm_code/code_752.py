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
    Detects if the synthesis route involves urea formation in the late stage (low depth).
    """
    urea_formation_detected = False
    max_depth_for_late_stage = 1  # Define what "late stage" means

    def dfs_traverse(node, depth=0):
        nonlocal urea_formation_detected

        if node["type"] == "reaction" and depth <= max_depth_for_late_stage:
            # Extract reaction SMILES
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains urea
                product_mol = Chem.MolFromSmiles(product_smiles)
                urea_pattern = Chem.MolFromSmarts("[#7]-[#6](=[#8])-[#7]")

                # Check if reactants contain isocyanate
                reactants = reactants_smiles.split(".")
                isocyanate_pattern = Chem.MolFromSmarts("[#8]=[#6]=[#7]")
                has_isocyanate = any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(isocyanate_pattern)
                    for r in reactants
                    if r
                )

                if product_mol and product_mol.HasSubstructMatch(urea_pattern) and has_isocyanate:
                    print(f"Late-stage urea formation detected at depth {depth}")
                    urea_formation_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return urea_formation_detected
