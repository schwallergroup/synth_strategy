#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


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
