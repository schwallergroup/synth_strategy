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
    This function detects the overall strategy: a linear synthesis featuring carboxylic acid
    protection-deprotection cycle with benzyl group, combined with sequential nitrogen functional
    group interconversions, and late-stage esterification via acid chloride activation.
    """
    print("Analyzing synthesis route for combined strategy...")

    # Check all individual strategies
    print("\nChecking for benzyl protection-deprotection strategy:")
    has_benzyl_protection = benzyl_protection_deprotection_strategy(route)
    print(f"Benzyl protection-deprotection strategy found: {has_benzyl_protection}")

    print("\nChecking for nitrogen functional group interconversions:")
    has_nitrogen_interconversion = nitrogen_functional_group_interconversion(route)
    print(f"Sequential nitrogen interconversions found: {has_nitrogen_interconversion}")

    print("\nChecking for acid chloride esterification:")
    has_acid_chloride_esterification = acid_chloride_activation_esterification(route)
    print(f"Late-stage acid chloride esterification found: {has_acid_chloride_esterification}")

    print("\nChecking for linear synthesis strategy:")
    is_linear = linear_synthesis_strategy(route)
    print(f"Linear synthesis strategy found: {is_linear}")

    # The combined strategy requires all individual strategies to be present
    combined_strategy_present = (
        has_benzyl_protection
        and has_nitrogen_interconversion
        and has_acid_chloride_esterification
        and is_linear
    )

    print("\nFinal result:")
    if combined_strategy_present:
        print(
            "Combined strategy detected: Linear synthesis with benzyl protection-deprotection, "
            "nitrogen functional group interconversions, and acid chloride esterification"
        )
    else:
        print("Combined strategy NOT detected - missing some required elements")

    return combined_strategy_present
