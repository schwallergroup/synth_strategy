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
    This function detects if the synthesis involves heterocycle formation as a key strategy.
    It combines the detection of heterocycle formation with the overall linear synthesis approach.
    """
    # Track synthesis stages and chemical transformations
    synthesis_data = {
        "heterocycle_formation": False,
        "heterocycle_formation_depth": float("inf"),
        "heterocycle_type": None,
        "max_depth": 0,
        "linear_structure": True,
        "oxime_reduction": False,
        "oxime_reduction_depth": float("inf"),
        "halide_to_alcohol": False,
        "halide_to_alcohol_depth": float("inf"),
        "other_fg_transformation": False,
        "other_fg_transformation_depth": float("inf"),
        "branching_factor": [],
    }

    # Traverse the synthesis route to analyze its structure and transformations
    dfs_traverse(route, 0, synthesis_data)

    # Check if heterocycle formation occurred in the middle of synthesis
    has_heterocycle_formation = heterocycle_formation_mid_synthesis(synthesis_data)

    # Check if synthesis follows a linear approach
    has_linear_strategy = linear_synthesis_strategy(synthesis_data)

    # Check for specific functional group transformations at appropriate stages
    has_fg_transformation = (
        oxime_reduction_late_stage(synthesis_data)
        or halide_to_alcohol_early_stage(synthesis_data)
        or other_fg_transformation_detected(synthesis_data)
    )

    # The strategy is present if we have heterocycle formation in a linear synthesis with FG transformations
    strategy_present = has_heterocycle_formation and has_linear_strategy and has_fg_transformation

    if strategy_present:
        print(
            f"Detected heterocycle formation strategy ({synthesis_data['heterocycle_type']}) with functional group transformations"
        )
    else:
        print(f"Missing conditions for heterocycle formation strategy:")
        print(f"  - Heterocycle formation: {has_heterocycle_formation}")
        print(f"  - Linear synthesis: {has_linear_strategy}")
        print(f"  - FG transformation: {has_fg_transformation}")

    return strategy_present
