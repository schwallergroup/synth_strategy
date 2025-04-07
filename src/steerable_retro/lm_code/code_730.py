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
    This function detects a protection-deprotection sequence involving acetate groups.
    """
    # Track protection and deprotection events
    protection_events = []
    deprotection_events = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants.split(".")]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol:
                    # Define SMARTS patterns
                    acetate_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][#6]")
                    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")

                    # Check for protection (alcohol → acetate)
                    alcohol_in_reactants = any(
                        r and r.HasSubstructMatch(alcohol_pattern) for r in reactants_mols if r
                    )
                    acetate_in_product = product_mol.HasSubstructMatch(acetate_pattern)

                    if alcohol_in_reactants and acetate_in_product:
                        protection_events.append(node.get("metadata", {}).get("ID", "unknown"))

                    # Check for deprotection (acetate → alcohol)
                    acetate_in_reactants = any(
                        r and r.HasSubstructMatch(acetate_pattern) for r in reactants_mols if r
                    )
                    alcohol_in_product = product_mol.HasSubstructMatch(alcohol_pattern)

                    if acetate_in_reactants and alcohol_in_product:
                        deprotection_events.append(node.get("metadata", {}).get("ID", "unknown"))

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Check if both protection and deprotection occurred
    has_protection_deprotection = len(protection_events) > 0 and len(deprotection_events) > 0

    if has_protection_deprotection:
        print(
            f"Detected protection-deprotection sequence: Protection at {protection_events}, Deprotection at {deprotection_events}"
        )

    return has_protection_deprotection
