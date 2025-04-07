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
    This function detects if the synthetic route involves benzyl ether formation
    as a protection strategy for phenols.
    """
    benzyl_ether_count = 0
    phenol_protection_events = 0

    def dfs_traverse(node):
        nonlocal benzyl_ether_count, phenol_protection_events

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for phenol to benzyl ether transformation
            phenol_pattern = Chem.MolFromSmarts("[c][OH]")
            benzyl_ether_pattern = Chem.MolFromSmarts("[c][O][CH2][c]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol:
                # Check if product has benzyl ether
                if product_mol.HasSubstructMatch(benzyl_ether_pattern):
                    benzyl_ether_count += 1

                    # Check if any reactant has phenol
                    for r_mol in reactant_mols:
                        if r_mol and r_mol.HasSubstructMatch(phenol_pattern):
                            phenol_protection_events += 1
                            print(f"Found phenol protection event: {rsmi}")
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # Strategy is present if we have at least one phenol protection event
    result = phenol_protection_events >= 1
    print(
        f"Benzyl ether protection strategy detected: {result} (protection events: {phenol_protection_events})"
    )
    return result
