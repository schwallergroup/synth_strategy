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
    This function detects a silyl protection-deprotection sequence in the synthesis.
    Specifically looking for O-Si bond formation followed by cleavage.
    """
    # Track protection and deprotection events
    protection_events = []  # Will store (depth, product_smiles, reactant_smiles)
    deprotection_events = []  # Will store (depth, product_smiles, reactant_smiles)

    def dfs_traverse(node, depth=0):
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Extract reactants and product
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]
                print(f"  Product: {product}")
                print(f"  Reactants: {reactants}")
            except Exception as e:
                print(f"Error extracting reactants/product: {e}")
                return

            # Check for silyl protection reaction - either by reaction type or by functional group changes
            if checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                print(f"Silyl protection found at depth {depth} (by reaction type)")
                protection_events.append((depth, product, reactants))
            # Alternative detection: check if alcohol in reactants and silyl group in product
            elif any(
                checker.check_fg("Alcohol", r) for r in reactants
            ) and checker.check_fg("Silyl protective group", product):
                # Verify this is likely a protection by checking for silyl chloride in reactants
                if any("Si" in r and "Cl" in r for r in reactants):
                    print(
                        f"Silyl protection found at depth {depth} (by functional group analysis)"
                    )
                    protection_events.append((depth, product, reactants))

            # Check for silyl deprotection reactions - either by reaction type or by functional group changes
            if (
                checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
                or checker.check_reaction(
                    "Alcohol deprotection from silyl ethers (double)", rsmi
                )
                or checker.check_reaction(
                    "Alcohol deprotection from silyl ethers (diol)", rsmi
                )
            ):
                print(f"Silyl deprotection found at depth {depth} (by reaction type)")
                deprotection_events.append((depth, product, reactants))
            # Alternative detection: check if silyl group in reactants and alcohol in product
            elif any(
                checker.check_fg("Silyl protective group", r) for r in reactants
            ) and checker.check_fg("Alcohol", product):
                print(
                    f"Silyl deprotection found at depth {depth} (by functional group analysis)"
                )
                deprotection_events.append((depth, product, reactants))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Protection events: {len(protection_events)}")
    print(f"Deprotection events: {len(deprotection_events)}")

    # Check if both protection and deprotection were found
    if not protection_events or not deprotection_events:
        print("Missing either protection or deprotection events")
        return False

    # In retrosynthesis, higher depth means earlier stage in forward synthesis
    # Protection should occur before deprotection in forward synthesis
    for prot_depth, prot_product, prot_reactants in protection_events:
        for deprot_depth, deprot_product, deprot_reactants in deprotection_events:
            # Check if protection occurs before deprotection in forward synthesis
            if prot_depth > deprot_depth:
                print(
                    f"Found potential sequence: protection at depth {prot_depth}, deprotection at depth {deprot_depth}"
                )

                # Check if the protected molecule is related to the deprotected one
                # In a proper sequence, the product of protection should be involved in deprotection
                for reactant in deprot_reactants:
                    if checker.check_fg(
                        "Silyl protective group", reactant
                    ) and checker.check_fg("Alcohol", deprot_product):
                        print("Silyl protection-deprotection strategy detected")
                        return True

    print("No valid silyl protection-deprotection sequence found")
    return False
