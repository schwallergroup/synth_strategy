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
    This function detects if the synthesis route involves Boc protection/deprotection sequence.
    """
    # Track protection and deprotection events with their depths
    protection_events = []
    deprotection_events = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Depth {depth}, Examining reaction: {rsmi}")

            # Extract reactants and product
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc protection reactions
            is_boc_protection = False
            if (
                checker.check_reaction("Boc amine protection", rsmi)
                or checker.check_reaction("Boc amine protection explicit", rsmi)
                or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                or checker.check_reaction("Boc amine protection of primary amine", rsmi)
            ):
                is_boc_protection = True
                print(f"Depth {depth}: Detected Boc protection reaction")

            # Alternative check for Boc protection by examining functional groups
            if not is_boc_protection:
                # Check if product contains Boc group and reactants contain amine
                has_boc_in_product = "C(=O)OC(C)(C)C" in product or checker.check_fg(
                    "Carbamic ester", product
                )
                has_amine_in_reactants = any(
                    checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                    for r in reactants
                )

                if has_boc_in_product and has_amine_in_reactants:
                    is_boc_protection = True
                    print(f"Depth {depth}: Detected Boc protection by functional group analysis")

            if is_boc_protection:
                protection_events.append((depth, rsmi))

            # Check for Boc deprotection reactions
            is_boc_deprotection = False
            if (
                checker.check_reaction("Boc amine deprotection", rsmi)
                or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
            ):
                is_boc_deprotection = True
                print(f"Depth {depth}: Detected Boc deprotection reaction")

            # Alternative check for Boc deprotection by examining functional groups
            if not is_boc_deprotection:
                # Check if reactants contain Boc group and product contains amine
                has_boc_in_reactants = any(
                    "C(=O)OC(C)(C)C" in r or checker.check_fg("Carbamic ester", r)
                    for r in reactants
                )
                has_amine_in_product = checker.check_fg(
                    "Primary amine", product
                ) or checker.check_fg("Secondary amine", product)

                if has_boc_in_reactants and has_amine_in_product:
                    is_boc_deprotection = True
                    print(f"Depth {depth}: Detected Boc deprotection by functional group analysis")

            if is_boc_deprotection:
                deprotection_events.append((depth, rsmi))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if both protection and deprotection occurred
    has_protection = len(protection_events) > 0
    has_deprotection = len(deprotection_events) > 0

    print(f"Boc protection events: {len(protection_events)}")
    for depth, rsmi in protection_events:
        print(f"  Depth {depth}: {rsmi}")

    print(f"Boc deprotection events: {len(deprotection_events)}")
    for depth, rsmi in deprotection_events:
        print(f"  Depth {depth}: {rsmi}")

    # In retrosynthetic analysis, deprotection should be encountered at a lower depth than protection
    # (since we're traversing backwards from the target)
    has_correct_sequence = False
    if has_protection and has_deprotection:
        # Get the minimum depth for each event type
        min_protection_depth = min([depth for depth, _ in protection_events])
        min_deprotection_depth = min([depth for depth, _ in deprotection_events])

        # In retrosynthesis, deprotection should be encountered before protection
        # (lower depth = later stage in forward synthesis)
        if min_deprotection_depth < min_protection_depth:
            has_correct_sequence = True
            print(
                f"Correct sequence: Deprotection (depth {min_deprotection_depth}) occurs before Protection (depth {min_protection_depth}) in retrosynthetic analysis"
            )
        else:
            print(
                f"Incorrect sequence: Protection (depth {min_protection_depth}) occurs before Deprotection (depth {min_deprotection_depth}) in retrosynthetic analysis"
            )

    result = has_protection and has_deprotection and has_correct_sequence
    print(f"Boc protection/deprotection sequence detected: {result}")
    return result
