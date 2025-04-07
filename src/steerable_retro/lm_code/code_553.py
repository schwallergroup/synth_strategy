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

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects a synthesis strategy involving alcohol protection with silyl groups.
    It checks for both protection and deprotection steps and verifies that the protected
    intermediate is used in at least one reaction before deprotection.
    """
    # Track protected molecules and their deprotection
    protected_molecules = {}  # product SMILES -> depth
    deprotected_molecules = {}  # product SMILES -> depth
    protection_reactions = []
    deprotection_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alcohol protection with silyl group
            is_protection = checker.check_reaction("Alcohol protection with silyl ethers", rsmi)
            print(f"Checking protection at depth {depth}: {rsmi}, is_protection={is_protection}")

            if is_protection:
                # Verify reactant has alcohol group
                for reactant in reactants:
                    has_alcohol = (
                        checker.check_fg("Primary alcohol", reactant)
                        or checker.check_fg("Secondary alcohol", reactant)
                        or checker.check_fg("Tertiary alcohol", reactant)
                        or checker.check_fg("Aromatic alcohol", reactant)
                    )

                    # Verify product has silyl group
                    has_silyl = checker.check_fg(
                        "Silyl protective group", product
                    ) or checker.check_fg("TMS ether protective group", product)

                    print(f"  Reactant {reactant[:30]}... has_alcohol={has_alcohol}")
                    print(f"  Product {product[:30]}... has_silyl={has_silyl}")

                    if has_alcohol and has_silyl:
                        protected_molecules[product] = depth
                        protection_reactions.append(rsmi)
                        print(f"Found alcohol protection with silyl group at depth {depth}")

            # Check for alcohol deprotection from silyl ethers
            is_deprotection = (
                checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
                or checker.check_reaction("Alcohol deprotection from silyl ethers (double)", rsmi)
                or checker.check_reaction("Alcohol deprotection from silyl ethers (diol)", rsmi)
            )

            print(
                f"Checking deprotection at depth {depth}: {rsmi}, is_deprotection={is_deprotection}"
            )

            if is_deprotection:
                # Verify reactant has silyl group
                for reactant in reactants:
                    has_silyl = checker.check_fg(
                        "Silyl protective group", reactant
                    ) or checker.check_fg("TMS ether protective group", reactant)

                    # Verify product has alcohol group
                    has_alcohol = (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                        or checker.check_fg("Aromatic alcohol", product)
                    )

                    print(f"  Reactant {reactant[:30]}... has_silyl={has_silyl}")
                    print(f"  Product {product[:30]}... has_alcohol={has_alcohol}")

                    if has_silyl and has_alcohol:
                        deprotected_molecules[product] = depth
                        deprotection_reactions.append(rsmi)
                        print(f"Found alcohol deprotection from silyl ethers at depth {depth}")

            # Fallback: Check for pattern-based protection/deprotection if reaction check failed
            if not is_protection and not is_deprotection:
                # Check if this might be a protection reaction not properly classified
                for reactant in reactants:
                    if (
                        checker.check_fg("Primary alcohol", reactant)
                        or checker.check_fg("Secondary alcohol", reactant)
                        or checker.check_fg("Tertiary alcohol", reactant)
                        or checker.check_fg("Aromatic alcohol", reactant)
                    ):

                        # Look for silyl reagents in other reactants
                        silyl_reagent_present = any("Si" in r for r in reactants if r != reactant)

                        if silyl_reagent_present and (
                            checker.check_fg("Silyl protective group", product)
                            or checker.check_fg("TMS ether protective group", product)
                        ):
                            print(
                                f"Found potential silyl protection via pattern matching at depth {depth}"
                            )
                            protected_molecules[product] = depth
                            protection_reactions.append(rsmi)

                # Check if this might be a deprotection reaction not properly classified
                for reactant in reactants:
                    if checker.check_fg("Silyl protective group", reactant) or checker.check_fg(
                        "TMS ether protective group", reactant
                    ):

                        if (
                            checker.check_fg("Primary alcohol", product)
                            or checker.check_fg("Secondary alcohol", product)
                            or checker.check_fg("Tertiary alcohol", product)
                            or checker.check_fg("Aromatic alcohol", product)
                        ):
                            print(
                                f"Found potential silyl deprotection via pattern matching at depth {depth}"
                            )
                            deprotected_molecules[product] = depth
                            deprotection_reactions.append(rsmi)

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Protected molecules: {len(protected_molecules)}")
    print(f"Deprotected molecules: {len(deprotected_molecules)}")
    print(f"Protection reactions: {len(protection_reactions)}")
    print(f"Deprotection reactions: {len(deprotection_reactions)}")

    # Check if we have both protection and deprotection steps
    if protection_reactions and deprotection_reactions:
        # Check if any protected molecule was used in intermediate steps
        # In retrosynthesis, protection (earlier) should have higher depth than deprotection (later)
        for prot_depth in protected_molecules.values():
            for deprot_depth in deprotected_molecules.values():
                print(f"Comparing depths: protection={prot_depth}, deprotection={deprot_depth}")
                if (
                    prot_depth > deprot_depth
                ):  # Protection happened before deprotection in retrosynthesis
                    print(
                        "Found complete alcohol protection-deprotection strategy with silyl ethers"
                    )
                    return True

    # If we have at least one protection step, consider it a partial strategy
    if protection_reactions:
        print(
            "Found alcohol protection with silyl ethers, but no clear protection-deprotection strategy"
        )
        return True

    # If we have at least one deprotection step, it's still a silyl strategy
    if deprotection_reactions:
        print("Found alcohol deprotection from silyl ethers, but no clear protection step")
        return True

    print("No alcohol protection/deprotection with silyl ethers found")
    return False
