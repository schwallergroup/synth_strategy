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
    This function detects reductive amination (aldehyde/ketone + amine → secondary/tertiary amine).
    """
    reductive_amination_found = False

    def dfs_traverse(node):
        nonlocal reductive_amination_found

        if node["type"] == "reaction" and not reductive_amination_found:
            # Extract reaction SMILES
            rsmi = node["metadata"]["rsmi"]

            # Check for reductive amination reactions directly
            if checker.check_reaction("Reductive amination with aldehyde", rsmi):
                print(f"Reductive amination with aldehyde detected: {rsmi}")
                reductive_amination_found = True
                return

            if checker.check_reaction("Reductive amination with ketone", rsmi):
                print(f"Reductive amination with ketone detected: {rsmi}")
                reductive_amination_found = True
                return

            if checker.check_reaction("Reductive amination with alcohol", rsmi):
                print(f"Reductive amination with alcohol detected: {rsmi}")
                reductive_amination_found = True
                return

            # If direct reaction check fails, try checking by functional groups
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for aldehyde/ketone and amine in reactants
            has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants_smiles)
            has_ketone = any(checker.check_fg("Ketone", r) for r in reactants_smiles)
            has_primary_amine = any(checker.check_fg("Primary amine", r) for r in reactants_smiles)
            has_secondary_amine = any(
                checker.check_fg("Secondary amine", r) for r in reactants_smiles
            )

            # Check for secondary/tertiary amine in product
            has_sec_amine_product = checker.check_fg("Secondary amine", product_smiles)
            has_tert_amine_product = checker.check_fg("Tertiary amine", product_smiles)

            # Verify reductive amination pattern
            if (
                (has_aldehyde or has_ketone)
                and (has_primary_amine or has_secondary_amine)
                and (has_sec_amine_product or has_tert_amine_product)
            ):
                print(f"Reductive amination detected by functional group analysis: {rsmi}")
                reductive_amination_found = True

        # Traverse children
        for child in node.get("children", []):
            if not reductive_amination_found:  # Stop traversal if already found
                dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return reductive_amination_found
