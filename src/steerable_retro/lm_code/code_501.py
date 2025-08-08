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
    Detects the complete strategy: cyclopropane formation in retrosynthesis,
    with heterocycle preservation and late-stage esterification.
    """
    # Check all component strategies
    has_cyclopropane_formation = cyclopropane_ring_formation_strategy(route)
    has_late_esterification = late_stage_esterification(route)
    has_heterocycle_preservation = heterocycle_preservation_strategy(route)

    print(f"Cyclopropane formation: {has_cyclopropane_formation}")
    print(f"Late-stage esterification: {has_late_esterification}")
    print(f"Heterocycle preservation: {has_heterocycle_preservation}")

    # The complete strategy requires all three components
    return has_cyclopropane_formation and has_late_esterification and has_heterocycle_preservation
