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
    Detects if the synthesis uses a masked nucleophile strategy where a nitrogen nucleophile
    is installed through sequential transformations
    """
    # This strategy combines the alcohol→bromide→azide→amine sequence with the presence
    # of a nitro group that's later reduced to an amine

    # Check for both component strategies
    has_leaving_group_sequence = amine_installation_via_leaving_group_sequence(route)
    has_nitro_reduction = nitro_reduction_final_step(route)

    # The strategy is present if both components are found
    if has_leaving_group_sequence and has_nitro_reduction:
        print("Found masked nucleophile strategy")
        return True
    return False
