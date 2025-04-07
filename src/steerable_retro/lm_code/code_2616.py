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
    Detects if the synthesis includes an early-stage sulfonylation (installation of sulfonyl group).
    """
    # Track if we found the pattern and at what depth
    found_sulfonylation = False
    sulfonylation_depth = -1
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_sulfonylation, sulfonylation_depth, max_depth

        # Track maximum depth to determine early vs late stage
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for sulfonylation pattern
                try:
                    product_mol = Chem.MolFromSmiles(product)

                    if product_mol:
                        # Sulfonyl group pattern
                        sulfonyl_pattern = Chem.MolFromSmarts("[#16](=[O])(=[O])")

                        if product_mol.HasSubstructMatch(sulfonyl_pattern):
                            # Check if reactants don't have the sulfonyl group
                            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                            if not any(
                                mol and mol.HasSubstructMatch(sulfonyl_pattern)
                                for mol in reactant_mols
                                if mol
                            ):
                                found_sulfonylation = True
                                sulfonylation_depth = depth
                                print(f"Found sulfonylation at depth {depth}")
                except:
                    print("Error processing SMILES in sulfonylation check")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if sulfonylation is found at a high depth (early stage)
    # We consider it early stage if it's in the last third of the synthesis depth
    return found_sulfonylation and sulfonylation_depth >= (max_depth * 2 / 3)
