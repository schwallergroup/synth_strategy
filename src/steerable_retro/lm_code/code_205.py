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
    This function detects the combined strategy: linear synthesis with multiple ether formations
    and protection/deprotection sequences, featuring late-stage incorporation of a cycloalkyl group
    via reductive amination.
    """
    # Check all individual strategies
    is_linear = linear_synthesis_strategy(route)
    has_multiple_ethers = multiple_ether_formation_strategy(route)
    has_protection_deprotection = protection_deprotection_strategy(route)
    has_late_cycloalkyl = late_stage_cycloalkyl_incorporation(route)
    has_reductive_amination = reductive_amination_strategy(route)

    # Combined strategy requires most of these elements
    strategy_count = sum(
        [
            is_linear,
            has_multiple_ethers,
            has_protection_deprotection,
            has_late_cycloalkyl,
            has_reductive_amination,
        ]
    )

    print(f"Strategy elements detected: {strategy_count}/5")
    print(f"- Linear synthesis: {is_linear}")
    print(f"- Multiple ether formations: {has_multiple_ethers}")
    print(f"- Protection/deprotection: {has_protection_deprotection}")
    print(f"- Late-stage cycloalkyl: {has_late_cycloalkyl}")
    print(f"- Reductive amination: {has_reductive_amination}")

    # Require at least 3 of the 5 elements to match the combined strategy
    return strategy_count >= 3
