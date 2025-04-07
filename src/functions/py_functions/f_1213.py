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
    This function detects the complete strategy: a linear synthesis with Boc protection-deprotection,
    sulfonamide formation, nitro reduction, amide formation, and late-stage N-alkylation to build
    a piperazine sulfonamide scaffold.
    """
    # Check all individual components
    has_boc_seq = has_boc_protection_deprotection_sequence(route)
    has_nitro_red = has_nitro_reduction_to_amine(route)
    has_sulfonamide = has_sulfonamide_formation(route)
    has_amide = has_amide_formation_from_acid_chloride(route)
    has_n_alkylation = has_late_stage_n_alkylation(route)
    has_linear = has_linear_synthesis_strategy(route)
    has_scaffold = has_piperazine_sulfonamide_scaffold(route)

    print(f"Boc protection-deprotection: {has_boc_seq}")
    print(f"Nitro reduction: {has_nitro_red}")
    print(f"Sulfonamide formation: {has_sulfonamide}")
    print(f"Amide formation: {has_amide}")
    print(f"Late-stage N-alkylation: {has_n_alkylation}")
    print(f"Linear synthesis: {has_linear}")
    print(f"Piperazine sulfonamide scaffold: {has_scaffold}")

    # The complete strategy requires the scaffold, linear synthesis, and at least one key reaction
    result = (
        has_scaffold
        and has_linear
        and (
            has_boc_seq
            or has_nitro_red
            or has_sulfonamide
            or has_amide
            or has_n_alkylation
        )
    )

    print(f"Complete piperazine sulfonamide strategy detected: {result}")
    return result
