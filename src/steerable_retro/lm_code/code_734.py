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
    This function detects the use of azide as an intermediate for amine synthesis.
    It looks for mesylate formation, azide substitution, and azide reduction sequence.
    """
    mesylate_formation = False
    azide_formation = False
    azide_reduction = False

    def dfs_traverse(node):
        nonlocal mesylate_formation, azide_formation, azide_reduction

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Checking reaction: {rsmi}")

            # Check for mesylate/tosylate/triflate formation (alcohol to leaving group)
            if any(
                checker.check_fg("Primary alcohol", r)
                or checker.check_fg("Secondary alcohol", r)
                or checker.check_fg("Tertiary alcohol", r)
                for r in reactants
            ) and (
                checker.check_fg("Mesylate", product)
                or checker.check_fg("Tosylate", product)
                or checker.check_fg("Triflate", product)
            ):
                print("Found mesylate/leaving group formation step")
                mesylate_formation = True

            # Check for azide substitution (mesylate/tosylate/triflate to azide)
            if any(
                checker.check_fg("Mesylate", r)
                or checker.check_fg("Tosylate", r)
                or checker.check_fg("Triflate", r)
                for r in reactants
            ) and checker.check_fg("Azide", product):
                print("Found azide formation step")
                azide_formation = True

            # Check for azide reduction (azide to amine) - any reduction method
            if any(checker.check_fg("Azide", r) for r in reactants) and checker.check_fg(
                "Primary amine", product
            ):
                print("Found azide reduction step")
                azide_reduction = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if all three steps are found
    result = mesylate_formation and azide_formation and azide_reduction
    print(f"Strategy detection result: {result}")
    print(
        f"Steps found - Mesylate formation: {mesylate_formation}, Azide formation: {azide_formation}, Azide reduction: {azide_reduction}"
    )

    return result
