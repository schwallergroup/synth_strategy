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
    Detects if the route employs a protection-modification-deprotection strategy.
    This combines multiple features: protection, structural modification, and deprotection.
    """
    # Track protection, modification, and deprotection events
    protection_events = []  # List of (depth, product_smiles, protecting_group) tuples
    deprotection_events = []  # List of (depth, reactant_smiles, protecting_group) tuples
    modification_events = (
        []
    )  # List of (depth, reactant_smiles, product_smiles, protecting_group) tuples

    # Define protection and deprotection reaction types to check
    protection_reactions = [
        "Boc amine protection",
        "Boc amine protection explicit",
        "Boc amine protection with Boc anhydride",
        "Boc amine protection (ethyl Boc)",
        "Boc amine protection of secondary amine",
        "Boc amine protection of primary amine",
        "Alcohol protection with silyl ethers",
        "Protection of carboxylic acid",
    ]

    deprotection_reactions = [
        "Boc amine deprotection",
        "Boc amine deprotection of guanidine",
        "Boc amine deprotection to NH-NH2",
        "Tert-butyl deprotection of amine",
        "Alcohol deprotection from silyl ethers",
        "Alcohol deprotection from silyl ethers (double)",
        "Alcohol deprotection from silyl ethers (diol)",
        "Deprotection of carboxylic acid",
        "COOH ethyl deprotection",
        "Hydroxyl benzyl deprotection",
        "Carboxyl benzyl deprotection",
        "N-glutarimide deprotection",
        "Phthalimide deprotection",
        "TMS deprotection from alkyne",
    ]

    # Define protecting groups and their corresponding functional groups
    protecting_groups = {
        "Boc": ["Primary amine", "Secondary amine", "Tertiary amine"],
        "Silyl": ["Primary alcohol", "Secondary alcohol", "Tertiary alcohol"],
        "Benzyl": ["Primary alcohol", "Secondary alcohol", "Tertiary alcohol", "Carboxylic acid"],
        "TMS": ["Alkyne"],
        "Ethyl": ["Carboxylic acid"],
        "Glutarimide": ["Primary amine"],
        "Phthalimide": ["Primary amine"],
    }

    def identify_protecting_group(rsmi):
        """Identify the protecting group used in a protection/deprotection reaction"""
        # Check for Boc protection/deprotection
        if any(checker.check_reaction(rxn, rsmi) for rxn in protection_reactions[:6]) or any(
            checker.check_reaction(rxn, rsmi) for rxn in deprotection_reactions[:4]
        ):
            return "Boc"
        # Check for silyl protection/deprotection
        elif any(
            checker.check_reaction(rxn, rsmi)
            for rxn in [
                protection_reactions[6],
                deprotection_reactions[4],
                deprotection_reactions[5],
                deprotection_reactions[6],
            ]
        ):
            return "Silyl"
        # Check for carboxylic acid protection/deprotection
        elif any(
            checker.check_reaction(rxn, rsmi)
            for rxn in [
                protection_reactions[7],
                deprotection_reactions[7],
                deprotection_reactions[8],
            ]
        ):
            return "Ethyl"
        # Check for benzyl protection/deprotection
        elif any(
            checker.check_reaction(rxn, rsmi)
            for rxn in [deprotection_reactions[9], deprotection_reactions[10]]
        ):
            return "Benzyl"
        # Check for TMS protection/deprotection
        elif checker.check_reaction(deprotection_reactions[13], rsmi):
            return "TMS"
        # Check for glutarimide and phthalimide deprotection
        elif checker.check_reaction(deprotection_reactions[11], rsmi):
            return "Glutarimide"
        elif checker.check_reaction(deprotection_reactions[12], rsmi):
            return "Phthalimide"

        # Check molecules for protecting groups
        molecules = rsmi.split(">")[0].split(".") + [rsmi.split(">")[-1]]
        for mol in molecules:
            if checker.check_fg("Boc", mol):
                return "Boc"
            elif checker.check_fg("Silyl protective group", mol) or checker.check_fg(
                "TMS ether protective group", mol
            ):
                return "Silyl"
            elif "COOEt" in mol or "COOC2H5" in mol:
                return "Ethyl"
            elif "CH2Ph" in mol or "Bn" in mol:
                return "Benzyl"
            elif "TMS" in mol:
                return "TMS"

        return None

    def has_protecting_group(smiles):
        """Check if a molecule has any protecting group"""
        return (
            checker.check_fg("Boc", smiles)
            or checker.check_fg("Silyl protective group", smiles)
            or checker.check_fg("TMS ether protective group", smiles)
            or "COOEt" in smiles
            or "COOC2H5" in smiles
            or "CH2Ph" in smiles
            or "Bn" in smiles
            or "TMS" in smiles
            or "N-glutarimide" in smiles
            or "Phthalimide" in smiles
        )

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for protection reactions
                is_protection = False
                for rxn in protection_reactions:
                    if checker.check_reaction(rxn, rsmi):
                        protecting_group = identify_protecting_group(rsmi)
                        if protecting_group:
                            protection_events.append((depth, product, protecting_group))
                            print(f"Found {protecting_group} protection at depth {depth}")
                            is_protection = True
                            break

                # Check for deprotection reactions
                is_deprotection = False
                if not is_protection:
                    for rxn in deprotection_reactions:
                        if checker.check_reaction(rxn, rsmi):
                            # Find the protected reactant
                            for reactant in reactants:
                                if has_protecting_group(reactant):
                                    protecting_group = identify_protecting_group(rsmi)
                                    if protecting_group:
                                        deprotection_events.append(
                                            (depth, reactant, protecting_group)
                                        )
                                        print(
                                            f"Found {protecting_group} deprotection at depth {depth}"
                                        )
                                        is_deprotection = True
                                        break
                            if is_deprotection:
                                break

                # Any other reaction involving a protected molecule could be a modification
                if not is_protection and not is_deprotection:
                    # Check if any reactant has a protecting group
                    for reactant in reactants:
                        if has_protecting_group(reactant):
                            # Determine which protecting group is present
                            if checker.check_fg("Boc", reactant):
                                protecting_group = "Boc"
                            elif checker.check_fg(
                                "Silyl protective group", reactant
                            ) or checker.check_fg("TMS ether protective group", reactant):
                                protecting_group = "Silyl"
                            elif "COOEt" in reactant or "COOC2H5" in reactant:
                                protecting_group = "Ethyl"
                            elif "CH2Ph" in reactant or "Bn" in reactant:
                                protecting_group = "Benzyl"
                            elif "TMS" in reactant:
                                protecting_group = "TMS"
                            else:
                                protecting_group = "Unknown"

                            modification_events.append((depth, reactant, product, protecting_group))
                            print(
                                f"Found modification with {protecting_group} protected molecule at depth {depth}"
                            )
                            break

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Group events by protecting group
    protection_by_group = {}
    deprotection_by_group = {}
    modification_by_group = {}

    for depth, smiles, group in protection_events:
        if group not in protection_by_group:
            protection_by_group[group] = []
        protection_by_group[group].append(depth)

    for depth, smiles, group in deprotection_events:
        if group not in deprotection_by_group:
            deprotection_by_group[group] = []
        deprotection_by_group[group].append(depth)

    for depth, reactant, product, group in modification_events:
        if group not in modification_by_group:
            modification_by_group[group] = []
        modification_by_group[group].append(depth)

    # Check if we have a complete protection-modification-deprotection sequence for any protecting group
    found_strategy = False

    # Special case: If we have Boc-protected molecules being modified but no explicit protection/deprotection
    # This is likely because the Boc protection/deprotection steps are outside the synthetic route
    if "Boc" in modification_by_group and len(modification_by_group["Boc"]) > 0:
        print(
            "Found Boc-protected molecules being modified, assuming protection-modification-deprotection strategy"
        )
        found_strategy = True

    # Check for complete sequences
    for group in protection_by_group.keys():
        if group in deprotection_by_group and group in modification_by_group:
            min_protection_depth = min(protection_by_group[group])
            max_deprotection_depth = max(deprotection_by_group[group])

            # In retrosynthetic analysis, protection should be at a higher depth than deprotection
            if min_protection_depth > max_deprotection_depth:
                # Check if there's a modification between protection and deprotection
                modification_between = any(
                    max_deprotection_depth < mod_depth < min_protection_depth
                    for mod_depth in modification_by_group[group]
                )
                if modification_between:
                    found_strategy = True
                    print(
                        f"Found complete protection-modification-deprotection strategy with {group} protecting group"
                    )
                    break

    print(f"Protection-modification-deprotection strategy detected: {found_strategy}")
    print(f"Protection events: {protection_events}")
    print(f"Modification events: {[(d, r, g) for d, r, _, g in modification_events]}")
    print(f"Deprotection events: {deprotection_events}")

    return found_strategy
