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


def main(data):
    """
    Check if halide to alcohol conversion occurred in the early stages of synthesis.
    Early stage is defined as in the first 25% of the maximum depth.
    """
    if not data["halide_to_alcohol"]:
        return False

    max_depth = data["max_depth"]
    halide_depth = data["halide_to_alcohol_depth"]

    # Define early stage as first 25% of max depth
    early_stage_threshold = max_depth * 0.25

    is_early_stage = halide_depth <= early_stage_threshold

    if is_early_stage:
        print(
            f"Halide to alcohol conversion occurred in early stage (depth {halide_depth} of {max_depth})"
        )
    else:
        print(
            f"Halide to alcohol conversion did NOT occur in early stage (depth {halide_depth} of {max_depth}, threshold: {early_stage_threshold})"
        )

    return is_early_stage
