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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    This function detects the overall strategy: convergent synthesis with tertiary alcohol formation,
    N-Boc protection, and incorporation of trifluoromethyl and pyridine groups.
    """
    # Check all individual strategies
    has_boc_protection = boc_protection_strategy(route)
    has_tertiary_alcohol = tertiary_alcohol_formation_strategy(route)
    has_trifluoromethyl = trifluoromethyl_containing_strategy(route)
    has_pyridine = pyridine_containing_strategy(route)
    is_convergent = convergent_synthesis_strategy(route)
    has_late_stage_alcohol = late_stage_alcohol_formation_strategy(route)

    print(f"Boc protection: {has_boc_protection}")
    print(f"Tertiary alcohol: {has_tertiary_alcohol}")
    print(f"Trifluoromethyl: {has_trifluoromethyl}")
    print(f"Pyridine: {has_pyridine}")
    print(f"Convergent: {is_convergent}")
    print(f"Late-stage alcohol: {has_late_stage_alcohol}")

    # Combined strategy requires core elements plus at least one from each category
    strategy_detected = (
        has_trifluoromethyl
        and has_pyridine
        and (has_tertiary_alcohol or has_late_stage_alcohol)
        and (has_boc_protection or is_convergent)
    )

    # Alternative detection if the main criteria aren't met but we have strong evidence
    if not strategy_detected and has_trifluoromethyl and has_pyridine and is_convergent:
        # Check if the target molecule has a tertiary alcohol
        if route["type"] == "mol" and checker.check_fg("Tertiary alcohol", route["smiles"]):
            strategy_detected = True
            print("Strategy detected based on target molecule analysis")

    print(f"Combined strategy detected: {strategy_detected}")
    return strategy_detected
