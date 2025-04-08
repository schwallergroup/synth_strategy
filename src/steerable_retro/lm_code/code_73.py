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
    This function detects if the synthetic route involves late-stage chlorination of a primary alcohol.
    """
    primary_alcohol_pattern = Chem.MolFromSmarts("[CH2][OH]")
    chloride_pattern = Chem.MolFromSmarts("[CH2][Cl]")
    late_stage_chlorination_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_chlorination_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if any reactant contains primary alcohol
            reactant_has_alcohol = False
            for r_smiles in reactants_smiles:
                try:
                    r_mol = Chem.MolFromSmiles(r_smiles)
                    if r_mol and r_mol.HasSubstructMatch(primary_alcohol_pattern):
                        reactant_has_alcohol = True
                        break
                except:
                    continue

            # Check if product contains chloride
            try:
                p_mol = Chem.MolFromSmiles(product_smiles)
                product_has_chloride = p_mol and p_mol.HasSubstructMatch(chloride_pattern)
            except:
                product_has_chloride = False

            # If reactant has primary alcohol and product has chloride, and it's at depth 0-1, it's a late-stage chlorination
            if reactant_has_alcohol and product_has_chloride and depth <= 1:
                print(f"Late-stage chlorination detected in reaction at depth {depth}: {rsmi}")
                late_stage_chlorination_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return late_stage_chlorination_detected
