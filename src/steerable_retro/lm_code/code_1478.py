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
    Detects the combined protection-activation-etherification strategy with late-stage sulfonamide formation.
    """
    # Check for all components of the strategy
    has_phthalimide = phthalimide_protection_strategy(route)
    has_tosylation = tosylation_activation_strategy(route)
    has_etherification = etherification_via_tosylate_strategy(route)
    has_late_sulfonamide = late_stage_sulfonamide_formation(route)

    # The full strategy requires at least 3 of the 4 components
    strategy_score = sum(
        [has_phthalimide, has_tosylation, has_etherification, has_late_sulfonamide]
    )
    full_strategy_detected = strategy_score >= 3

    print(
        f"Protection-activation-etherification strategy detected: {full_strategy_detected} (score: {strategy_score}/4)"
    )
    return full_strategy_detected
