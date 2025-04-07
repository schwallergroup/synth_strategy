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
    This function detects the combined strategy of morpholine ring manipulation,
    SNAr reactions, and Boc protection/deprotection.
    """
    print("Checking for morpholine ring manipulation...")
    morpholine_result = has_morpholine_ring_manipulation(route)
    print(f"Morpholine ring manipulation found: {morpholine_result}")

    print("Checking for multiple SNAr reactions...")
    snar_result = has_multiple_snar_reactions(route)
    print(f"Multiple SNAr reactions found: {snar_result}")

    print("Checking for Boc protection/deprotection...")
    boc_result = has_boc_protection_deprotection(route)
    print(f"Boc protection/deprotection found: {boc_result}")

    result = morpholine_result and snar_result and boc_result
    print(f"Overall strategy detection result: {result}")

    return result
