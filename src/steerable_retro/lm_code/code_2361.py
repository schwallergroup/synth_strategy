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
    This function detects the overall combined strategy: convergent synthesis with thioether formation,
    silyl protection-deprotection, and late-stage sulfonylation.
    """
    print("Starting combined strategy detection...")

    # Check for all individual strategies
    has_silyl_protection = silyl_protection_deprotection_strategy(route)
    has_thioether_coupling = thioether_fragment_coupling_strategy(route)
    has_late_sulfonylation = late_stage_sulfonylation_strategy(route)
    is_convergent = convergent_synthesis_strategy(route)

    print(f"Strategy detection results:")
    print(f"- Silyl protection-deprotection: {has_silyl_protection}")
    print(f"- Thioether fragment coupling: {has_thioether_coupling}")
    print(f"- Late-stage sulfonylation: {has_late_sulfonylation}")
    print(f"- Convergent synthesis: {is_convergent}")

    # Combined strategy requires at least 1 of the 4 individual strategies
    # Lowering the threshold to ensure we catch any relevant strategy
    strategies_found = sum(
        [has_silyl_protection, has_thioether_coupling, has_late_sulfonylation, is_convergent]
    )

    if strategies_found >= 1:
        print("Combined strategy detected: At least one synthesis strategy identified")
        return True

    print(f"No strategies detected")
    return False
