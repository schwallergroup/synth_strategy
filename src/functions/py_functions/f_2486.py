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
    Detects if the synthesis uses a late-stage amide formation strategy
    where an amide bond is formed in the final or penultimate step.
    """
    amide_formation_detected = False
    late_stage_depth = 1  # Consider depth 0 or 1 as late stage

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_detected

        if node["type"] == "reaction" and depth <= late_stage_depth:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for amide formation
            product_mol = Chem.MolFromSmiles(product_smiles)
            amide_pattern = Chem.MolFromSmarts("[C;$(C=O)][N;!$(N=*)]")

            # Check if product has amide but at least one reactant doesn't
            if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                # Count amide bonds in product
                amide_matches_product = len(
                    product_mol.GetSubstructMatches(amide_pattern)
                )

                # Count amide bonds in reactants
                total_amide_matches_reactants = 0
                for r_smiles in reactants_smiles:
                    r_mol = Chem.MolFromSmiles(r_smiles)
                    if r_mol:
                        total_amide_matches_reactants += len(
                            r_mol.GetSubstructMatches(amide_pattern)
                        )

                # If product has more amide bonds than reactants combined, amide formation occurred
                if amide_matches_product > total_amide_matches_reactants:
                    amide_formation_detected = True
                    print(f"Amide formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return amide_formation_detected
