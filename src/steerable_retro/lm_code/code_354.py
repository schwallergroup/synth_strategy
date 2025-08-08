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
    This function detects if the route includes a reductive amination step.
    """
    reductive_amination_found = False

    def dfs_traverse(node, depth=0):
        nonlocal reductive_amination_found

        # Debug information about current node
        indent = "  " * depth
        if node["type"] == "mol":
            print(f"{indent}Examining molecule node: {node.get('smiles', 'No SMILES')[:30]}...")
        else:
            print(f"{indent}Examining reaction node")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            print(f"{indent}Analyzing reaction: {rsmi[:50]}...")

            # Check for reductive amination reactions using the checker function
            if (
                checker.check_reaction("reductive amination with aldehyde", rsmi)
                or checker.check_reaction("reductive amination with ketone", rsmi)
                or checker.check_reaction("reductive amination with alcohol", rsmi)
                or checker.check_reaction("Mignonac reaction", rsmi)
            ):

                print(
                    f"{indent}✓ Detected reductive amination step via direct reaction check: {rsmi}"
                )
                reductive_amination_found = True
                return

            # Alternative check if the specific reaction types are not detected
            if not reductive_amination_found:
                try:
                    # Extract reactants and product
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"{indent}Checking functional groups in reactants and products")

                    # Check for carbonyl compounds in reactants
                    has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants)
                    has_ketone = any(checker.check_fg("Ketone", r) for r in reactants)
                    has_formaldehyde = any(checker.check_fg("Formaldehyde", r) for r in reactants)

                    # Check for amines in reactants
                    has_primary_amine = any(checker.check_fg("Primary amine", r) for r in reactants)
                    has_secondary_amine = any(
                        checker.check_fg("Secondary amine", r) for r in reactants
                    )
                    has_ammonia = any(
                        "N" in r and len(r) <= 3 for r in reactants
                    )  # Simple check for ammonia

                    # Check for amine products
                    has_secondary_amine_product = checker.check_fg("Secondary amine", product)
                    has_tertiary_amine_product = checker.check_fg("Tertiary amine", product)

                    # Print debug info about functional groups
                    print(
                        f"{indent}Reactants - Aldehyde: {has_aldehyde}, Ketone: {has_ketone}, Formaldehyde: {has_formaldehyde}"
                    )
                    print(
                        f"{indent}Reactants - Primary amine: {has_primary_amine}, Secondary amine: {has_secondary_amine}, Ammonia: {has_ammonia}"
                    )
                    print(
                        f"{indent}Product - Secondary amine: {has_secondary_amine_product}, Tertiary amine: {has_tertiary_amine_product}"
                    )

                    # If we have the right reactants and products, it's likely a reductive amination
                    if (
                        (has_aldehyde or has_ketone or has_formaldehyde)
                        and (has_primary_amine or has_secondary_amine or has_ammonia)
                        and (has_secondary_amine_product or has_tertiary_amine_product)
                    ):

                        print(
                            f"{indent}✓ Detected reductive amination based on functional groups: {rsmi}"
                        )
                        reductive_amination_found = True
                        return

                    # Check for imine intermediates that might indicate reductive amination
                    if not reductive_amination_found:
                        has_imine_reactant = any(
                            checker.check_fg("Substituted imine", r)
                            or checker.check_fg("Unsubstituted imine", r)
                            for r in reactants
                        )

                        if has_imine_reactant and (
                            has_secondary_amine_product or has_tertiary_amine_product
                        ):
                            print(
                                f"{indent}✓ Detected reductive amination via imine reduction: {rsmi}"
                            )
                            reductive_amination_found = True
                            return

                except Exception as e:
                    print(f"{indent}Error analyzing reaction: {e}")
                    print(f"{indent}Reaction SMILES: {rsmi}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting route traversal to find reductive amination steps...")
    dfs_traverse(route)

    print(f"Reductive amination found in route: {reductive_amination_found}")
    return reductive_amination_found
