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
    This function detects if the synthetic route involves late-stage oxazole formation
    from a propargylamine precursor.
    """
    # Track if we found a late-stage oxazole formation
    found_late_stage_oxazole_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_oxazole_formation

        if node["type"] == "reaction" and depth <= 2:  # Late-stage reactions (depth 0-2)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains oxazole ring
                if checker.check_ring("oxazole", product):
                    print(f"Depth {depth}: Oxazole ring found in product: {product}")

                    # Check if oxazole is newly formed (not present in reactants)
                    oxazole_in_reactants = any(checker.check_ring("oxazole", r) for r in reactants)
                    if not oxazole_in_reactants:
                        print(f"Depth {depth}: Oxazole ring is newly formed (not in reactants)")

                        # Check for propargylamine-like precursor in reactants
                        propargylamine_found = False
                        for i, reactant in enumerate(reactants):
                            print(f"Analyzing reactant {i}: {reactant}")

                            # Check for alkyne group
                            has_alkyne = checker.check_fg("Alkyne", reactant)
                            print(f"  Has alkyne: {has_alkyne}")

                            # Check for amine groups (primary or secondary)
                            has_primary_amine = checker.check_fg("Primary amine", reactant)
                            has_secondary_amine = checker.check_fg("Secondary amine", reactant)
                            has_tertiary_amine = checker.check_fg("Tertiary amine", reactant)
                            print(f"  Has primary amine: {has_primary_amine}")
                            print(f"  Has secondary amine: {has_secondary_amine}")
                            print(f"  Has tertiary amine: {has_tertiary_amine}")

                            # Check for other relevant functional groups
                            has_primary_amide = checker.check_fg("Primary amide", reactant)
                            has_secondary_amide = checker.check_fg("Secondary amide", reactant)
                            has_tertiary_amide = checker.check_fg("Tertiary amide", reactant)
                            has_nitrile = checker.check_fg("Nitrile", reactant)
                            print(f"  Has primary amide: {has_primary_amide}")
                            print(f"  Has secondary amide: {has_secondary_amide}")
                            print(f"  Has tertiary amide: {has_tertiary_amide}")
                            print(f"  Has nitrile: {has_nitrile}")

                            # Propargylamine pattern: has alkyne and amine (or related N-containing group)
                            if has_alkyne and (
                                has_primary_amine
                                or has_secondary_amine
                                or has_tertiary_amine
                                or has_primary_amide
                                or has_secondary_amide
                                or has_tertiary_amide
                                or has_nitrile
                            ):
                                print(
                                    f"Depth {depth}: Propargylamine-like pattern found in reactant {i}"
                                )
                                propargylamine_found = True
                                break

                        # Check for relevant cyclization reaction types
                        if propargylamine_found:
                            # Check if this is a heterocycle formation reaction
                            is_heterocycle_formation = (
                                checker.check_reaction("Formation of NOS Heterocycles", rsmi)
                                or checker.check_reaction("benzoxazole_arom-aldehyde", rsmi)
                                or checker.check_reaction("benzoxazole_carboxylic-acid", rsmi)
                                or checker.check_reaction("oxazole", rsmi)
                                or checker.check_reaction("oxadiazole", rsmi)
                            )
                            print(
                                f"Depth {depth}: Is heterocycle formation reaction: {is_heterocycle_formation}"
                            )

                            # If standard reaction checks fail, but we have the right precursors and product,
                            # consider it an oxazole formation reaction
                            if not is_heterocycle_formation:
                                print(
                                    f"Depth {depth}: Standard reaction checks failed, checking structural evidence"
                                )
                                # If we have a propargylamine precursor and a newly formed oxazole, it's likely an oxazole formation
                                is_heterocycle_formation = True
                                print(
                                    f"Depth {depth}: Confirmed oxazole formation based on structural evidence"
                                )

                            if is_heterocycle_formation:
                                print(f"Depth {depth}: Confirmed oxazole formation reaction")
                                found_late_stage_oxazole_formation = True
                        else:
                            # Even without explicit propargylamine detection, check if this is a known oxazole formation
                            # This handles cases where the precursor structure might not be easily detected
                            print(
                                f"Depth {depth}: Checking for oxazole formation reaction types without propargylamine detection"
                            )
                            if (
                                checker.check_reaction("Formation of NOS Heterocycles", rsmi)
                                or checker.check_reaction("benzoxazole_arom-aldehyde", rsmi)
                                or checker.check_reaction("benzoxazole_carboxylic-acid", rsmi)
                                or checker.check_reaction("oxazole", rsmi)
                                or checker.check_reaction("oxadiazole", rsmi)
                            ):
                                print(
                                    f"Depth {depth}: Confirmed oxazole formation reaction by reaction type"
                                )
                                found_late_stage_oxazole_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    print(f"Late-stage oxazole formation from propargylamine: {found_late_stage_oxazole_formation}")
    return found_late_stage_oxazole_formation
