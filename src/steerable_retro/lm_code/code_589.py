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
    Check if heterocycle formation occurred in the middle stages of synthesis.
    Middle stage is defined as between 25% and 75% of the maximum depth.
    """
    if not data["heterocycle_formation"]:
        return False

    max_depth = data["max_depth"]
    heterocycle_depth = data["heterocycle_formation_depth"]

    # Define middle stage as between 25% and 75% of max depth
    lower_bound = max_depth * 0.25
    upper_bound = max_depth * 0.75

    is_mid_synthesis = lower_bound <= heterocycle_depth <= upper_bound

    if is_mid_synthesis:
        print(
            f"Heterocycle formation occurred in middle stage (depth {heterocycle_depth} of {max_depth})"
        )
    else:
        print(
            f"Heterocycle formation did NOT occur in middle stage (depth {heterocycle_depth} of {max_depth}, bounds: {lower_bound}-{upper_bound})"
        )

    return is_mid_synthesis
