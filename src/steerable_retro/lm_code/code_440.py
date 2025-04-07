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
    This function detects the combined strategy of multiple amide formations,
    thiazole ring formation, nitro reduction, and nitrile preservation in a linear synthesis.
    """
    # Check all individual strategies
    has_multiple_amides = multiple_amide_formations_strategy(route)
    has_thiazole_formation = thiazole_ring_formation_strategy(route)
    has_nitro_reduction = nitro_reduction_strategy(route)
    preserves_nitrile = nitrile_preservation_strategy(route)
    has_trifluoromethyl = trifluoromethyl_containing_strategy(route)
    is_linear = linear_synthesis_strategy(route)

    print(f"Multiple amide formations: {has_multiple_amides}")
    print(f"Thiazole ring formation: {has_thiazole_formation}")
    print(f"Nitro reduction: {has_nitro_reduction}")
    print(f"Nitrile preservation: {preserves_nitrile}")
    print(f"Trifluoromethyl containing: {has_trifluoromethyl}")
    print(f"Linear synthesis: {is_linear}")

    # Combined strategy requires all individual strategies to be present
    result = (
        has_multiple_amides
        and has_thiazole_formation
        and has_nitro_reduction
        and preserves_nitrile
        and has_trifluoromethyl
        and is_linear
    )

    print(f"Combined strategy detected: {result}")
    return result
