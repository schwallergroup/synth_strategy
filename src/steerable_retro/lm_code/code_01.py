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
    Detects a synthetic strategy involving biaryl formation via Suzuki coupling.
    """
    has_suzuki_coupling = False

    def dfs_traverse(node):
        nonlocal has_suzuki_coupling

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants = parts[0].split(".")
            product = parts[2]

            # Check for Suzuki coupling
            boronic_acid_pattern = re.compile(r"[cB]|B\(O\)|OB\(O\)")
            halogen_pattern = re.compile(r"[cBrI]|Br|I|Cl")

            has_boronic_acid = any(boronic_acid_pattern.search(r) for r in reactants)
            has_halogen = any(halogen_pattern.search(r) for r in reactants)

            if has_boronic_acid and has_halogen:
                # Further verify it's likely a Suzuki by checking for aromatic rings
                if any("c" in r for r in reactants):
                    has_suzuki_coupling = True
                    print(f"Detected Suzuki coupling: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Biaryl formation via Suzuki: {has_suzuki_coupling}")

    return has_suzuki_coupling
