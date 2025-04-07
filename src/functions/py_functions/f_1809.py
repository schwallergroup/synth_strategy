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
    Detects the complete strategy: convergent synthesis with piperazine scaffold,
    Boc protection/deprotection, alcohol activation, and late-stage amide formation.
    """
    # Check all individual strategies
    convergent = convergent_synthesis_with_piperazine(route)
    boc_sequence = boc_protection_deprotection_sequence(route)
    alcohol_activation = alcohol_to_bromide_activation(route)
    late_amide = late_stage_amide_formation(route)

    # The full strategy requires all components
    full_strategy_present = (
        convergent and boc_sequence and alcohol_activation and late_amide
    )

    if full_strategy_present:
        print("Found complete convergent piperazine strategy with all components")

    return full_strategy_present
