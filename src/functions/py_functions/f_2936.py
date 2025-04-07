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
    This function detects the combined strategy of convergent synthesis with
    late-stage sulfonamide formation.
    """
    # Check for convergent synthesis
    is_convergent = check_convergent_synthesis(route)

    # Check for late-stage sulfonamide formation
    has_late_sulfonamide = check_late_stage_sulfonamide(route)

    # The combined strategy requires both conditions to be true
    print(
        f"Is convergent: {is_convergent}, Has late sulfonamide: {has_late_sulfonamide}"
    )
    return is_convergent and has_late_sulfonamide
