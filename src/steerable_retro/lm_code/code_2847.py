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
    Detects the complete strategy: convergent synthesis with piperazine scaffold,
    Boc protection/deprotection, alcohol activation, and late-stage amide formation.
    """
    # Check all individual strategies
    convergent = convergent_synthesis_with_piperazine(route)
    boc_sequence = boc_protection_deprotection_sequence(route)
    alcohol_activation = alcohol_to_bromide_activation(route)
    late_amide = late_stage_amide_formation(route)

    # The full strategy requires all components
    full_strategy_present = convergent and boc_sequence and alcohol_activation and late_amide

    if full_strategy_present:
        print("Found complete convergent piperazine strategy with all components")

    return full_strategy_present
