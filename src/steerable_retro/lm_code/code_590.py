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
    Check if the synthesis follows a predominantly linear approach.
    A linear synthesis typically has a low average branching factor.
    """
    if not data["branching_factor"]:
        return True  # Default to True if no branching data

    # Calculate average branching factor
    avg_branching = sum(data["branching_factor"]) / len(data["branching_factor"])

    # Consider synthesis linear if average branching factor is close to 1
    is_linear = avg_branching <= 1.5

    if is_linear:
        print(f"Synthesis follows a linear approach (avg branching: {avg_branching:.2f})")
    else:
        print(f"Synthesis does NOT follow a linear approach (avg branching: {avg_branching:.2f})")

    return is_linear
