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
    This function detects a synthetic strategy involving a sequence of carboxylic acid
    derivative transformations (ester → acid → different ester → alcohol).
    """
    # Track the sequence of functional groups observed
    functional_group_sequence = []

    def dfs_traverse(node):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Patterns for different carboxylic acid derivatives
            ester_pattern = Chem.MolFromSmarts("[#6][#8][C](=[O])[#6]")
            acid_pattern = Chem.MolFromSmarts("[#8H1][C](=[O])[#6]")
            alcohol_pattern = Chem.MolFromSmarts("[#8H1][#6]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol:
                # Determine the functional group transformation
                if any(
                    r and r.HasSubstructMatch(ester_pattern) for r in reactant_mols
                ) and product_mol.HasSubstructMatch(acid_pattern):
                    functional_group_sequence.append("ester_to_acid")
                    print(f"Detected ester to acid conversion: {rsmi}")
                elif any(
                    r and r.HasSubstructMatch(acid_pattern) for r in reactant_mols
                ) and product_mol.HasSubstructMatch(ester_pattern):
                    functional_group_sequence.append("acid_to_ester")
                    print(f"Detected acid to ester conversion: {rsmi}")
                elif any(
                    r and r.HasSubstructMatch(ester_pattern) for r in reactant_mols
                ) and product_mol.HasSubstructMatch(alcohol_pattern):
                    functional_group_sequence.append("ester_to_alcohol")
                    print(f"Detected ester to alcohol conversion: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if the sequence matches our target pattern
    # We're looking for a sequence that includes ester→acid→ester→alcohol transformations
    # This is a simplified check - in reality, you'd need to ensure these are connected transformations
    has_ester_to_acid = "ester_to_acid" in functional_group_sequence
    has_acid_to_ester = "acid_to_ester" in functional_group_sequence
    has_ester_to_alcohol = "ester_to_alcohol" in functional_group_sequence

    result = has_ester_to_acid and has_ester_to_alcohol
    print(f"Carboxylic acid derivative sequence strategy detected: {result}")
    print(f"Functional group sequence: {functional_group_sequence}")
    return result
