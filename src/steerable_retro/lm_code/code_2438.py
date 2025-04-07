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

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects a Boc protection/deprotection sequence in the synthesis.
    Returns True if a Boc group is added early in the synthesis and removed later.
    """
    # Track protection and deprotection events with depth information
    protection_events = []  # Will store (depth, molecule_smiles, atom_indices)
    deprotection_events = []  # Will store (depth, molecule_smiles, atom_indices)

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc protection using checker function
            if (
                checker.check_reaction("Boc amine protection", rsmi)
                or checker.check_reaction("Boc amine protection explicit", rsmi)
                or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                or checker.check_reaction("Boc amine protection of primary amine", rsmi)
            ):
                print(f"Found Boc protection reaction at depth {depth}: {rsmi}")
                # Store the protection event with depth and product molecule
                protection_events.append((depth, product))

            # Check for Boc deprotection using checker function
            elif (
                checker.check_reaction("Boc amine deprotection", rsmi)
                or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
            ):
                print(f"Found Boc deprotection reaction at depth {depth}: {rsmi}")
                # Store the deprotection event with depth and reactant molecule
                # In deprotection, the Boc group is on the reactant
                for reactant in reactants:
                    if checker.check_fg("Carbamic ester", reactant):
                        deprotection_events.append((depth, reactant))
                        break

            # Fallback to pattern matching if checker doesn't identify the reactions
            else:
                # Check for Boc protection patterns
                boc_reagents = ["C(C)(C)OC(=O)", "OC(=O)OC(C)(C)C", "(Boc)2O", "Boc2O"]
                if any(
                    any(boc_r in r for boc_r in boc_reagents) for r in reactants
                ) and checker.check_fg("Carbamic ester", product):
                    # Verify an amine in reactants is being protected
                    if any(
                        checker.check_fg(fg, r)
                        for r in reactants
                        for fg in ["Primary amine", "Secondary amine"]
                    ):
                        print(
                            f"Found Boc protection reaction (pattern match) at depth {depth}: {rsmi}"
                        )
                        protection_events.append((depth, product))

                # Check for Boc deprotection patterns
                elif any(checker.check_fg("Carbamic ester", r) for r in reactants):
                    # Verify an amine in the product (deprotected)
                    if any(
                        checker.check_fg(fg, product)
                        for fg in ["Primary amine", "Secondary amine", "Tertiary amine"]
                    ):
                        print(
                            f"Found Boc deprotection reaction (pattern match) at depth {depth}: {rsmi}"
                        )
                        for reactant in reactants:
                            if checker.check_fg("Carbamic ester", reactant):
                                deprotection_events.append((depth, reactant))
                                break

        # Traverse children (going backward in synthesis)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have both protection and deprotection events
    has_protection = len(protection_events) > 0
    has_deprotection = len(deprotection_events) > 0

    print(f"Boc protection events: {len(protection_events)}")
    print(f"Boc deprotection events: {len(deprotection_events)}")

    # Check if protection happens at a higher depth (earlier in synthesis)
    # than deprotection (later in synthesis)
    valid_sequence = False
    if has_protection and has_deprotection:
        # Get the minimum depth of protection (earliest protection)
        min_protection_depth = min(depth for depth, _ in protection_events)
        # Get the maximum depth of deprotection (latest deprotection)
        max_deprotection_depth = max(depth for depth, _ in deprotection_events)

        # In a valid sequence, protection should happen earlier (higher depth)
        # than deprotection (lower depth)
        valid_sequence = min_protection_depth > max_deprotection_depth
        print(f"Min protection depth: {min_protection_depth}")
        print(f"Max deprotection depth: {max_deprotection_depth}")
        print(f"Valid sequence: {valid_sequence}")

    # Return True if we have both protection and deprotection in the correct sequence
    return has_protection and has_deprotection and valid_sequence
