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
    Check if oxime reduction occurred in the late stages of synthesis.
    Late stage is defined as in the last 25% of the maximum depth.
    """
    if not data["oxime_reduction"]:
        return False

    max_depth = data["max_depth"]
    oxime_depth = data["oxime_reduction_depth"]

    # Define late stage as last 25% of max depth
    late_stage_threshold = max_depth * 0.75

    is_late_stage = oxime_depth >= late_stage_threshold

    if is_late_stage:
        print(f"Oxime reduction occurred in late stage (depth {oxime_depth} of {max_depth})")
    else:
        print(
            f"Oxime reduction did NOT occur in late stage (depth {oxime_depth} of {max_depth}, threshold: {late_stage_threshold})"
        )

    return is_late_stage
