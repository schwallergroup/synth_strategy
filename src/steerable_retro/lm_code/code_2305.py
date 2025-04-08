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
    Detects the overall strategy: linear synthesis with late-stage N-alkylation,
    featuring a protection-deprotection sequence and ester interconversion.
    """
    # Check for all the individual components of the strategy
    has_protection = has_protection_deprotection_sequence(route)
    has_ester_changes = has_ester_interconversion_sequence(route)
    has_n_alkylation = has_late_stage_n_alkylation(route)
    is_linear = check_linear_synthesis(route)

    print(
        f"Strategy components: Protection-deprotection: {has_protection}, Ester interconversion: {has_ester_changes}, Late-stage N-alkylation: {has_n_alkylation}, Linear synthesis: {is_linear}"
    )

    # The strategy is present if all key elements are detected
    if has_protection and has_ester_changes and has_n_alkylation and is_linear:
        print(
            "Complete strategy detected: Linear synthesis with late-stage N-alkylation, protection-deprotection sequence, and ester interconversion"
        )
        return True

    # If not all components are present, check if we have at least 3 out of 4
    if sum([has_protection, has_ester_changes, has_n_alkylation, is_linear]) >= 3:
        print("Partial strategy detected: At least 3 out of 4 components found")
        return True

    return False
