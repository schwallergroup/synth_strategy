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
    This function detects the overall strategy: sequential Sonogashira couplings to construct
    a diarylacetylene scaffold, followed by nitro reduction and late-stage cyclization.
    """
    # Check for diarylacetylene scaffold formation via Sonogashira
    has_diarylacetylene, diarylacetylene_product = diarylacetylene_scaffold_strategy(route)

    if not has_diarylacetylene:
        print("Strategy not found: Missing diarylacetylene scaffold formation")
        return False

    # Check for nitro reduction to amine, related to the diarylacetylene scaffold
    has_nitro_reduction, amine_product = nitro_reduction_to_amine_strategy(
        route, diarylacetylene_product
    )

    if not has_nitro_reduction:
        print("Strategy not found: Missing nitro reduction related to diarylacetylene scaffold")
        return False

    # Check for late-stage heterocyclization involving the amine from nitro reduction
    has_late_cyclization = late_stage_heterocyclization_strategy(route, amine_product)

    if not has_late_cyclization:
        print("Strategy not found: Missing late-stage heterocyclization involving the amine")
        return False

    # The overall strategy requires all three components in the correct sequence
    print(
        "Detected overall strategy: Linear diarylacetylene construction with nitro reduction and late heterocyclization"
    )
    return True
