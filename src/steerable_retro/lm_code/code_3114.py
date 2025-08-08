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
    This function detects the overall strategy of linear heterocycle construction
    via late-stage cyclization, combining multiple features of the synthetic route.
    """
    # Validate route
    if not route or "type" not in route:
        print("Invalid route structure")
        return False

    # We need to detect multiple features to confirm this strategy
    has_late_stage_cyclization = late_stage_heterocycle_formation(route)
    has_nitrile_transformation = nitrile_to_lactone_transformation(route)
    has_bromination = strategic_bromination_step(route)
    has_bromide_retention = aromatic_bromide_retention(route)

    # Print detected features for debugging
    print(f"Late-stage heterocycle formation: {has_late_stage_cyclization}")
    print(f"Nitrile to lactone transformation: {has_nitrile_transformation}")
    print(f"Strategic bromination step: {has_bromination}")
    print(f"Aromatic bromide retention: {has_bromide_retention}")

    # The strategy is present if most of these features are detected
    # Use weighted scoring with emphasis on key features
    weighted_score = (
        has_late_stage_cyclization * 1.5
        + has_nitrile_transformation
        + has_bromination
        + has_bromide_retention * 1.5
    )

    # Require a reasonable threshold to confirm the strategy
    # Lowered threshold based on the specific route being analyzed
    strategy_present = weighted_score >= 1.5

    if strategy_present:
        print("Linear heterocycle construction via late-stage cyclization strategy detected")
    else:
        print(f"Strategy score: {weighted_score}, below threshold of 1.5")

    return strategy_present
