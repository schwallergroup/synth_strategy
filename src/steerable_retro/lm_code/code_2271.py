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
    This function detects a protection-deprotection sequence in the synthesis,
    specifically looking for Boc and Cbz protection groups.
    """
    # Track protection/deprotection events with molecule information
    protected_molecules = {
        "boc": {},  # Will store {molecule_smiles: depth}
        "cbz": {},  # Will store {molecule_smiles: depth}
    }
    deprotected_molecules = {
        "boc": {},  # Will store {molecule_smiles: depth}
        "cbz": {},  # Will store {molecule_smiles: depth}
    }

    # Add depth information to the route
    def add_depth(node, current_depth=0):
        node["depth"] = current_depth
        for child in node.get("children", []):
            add_depth(child, current_depth + 1)

    add_depth(route)

    def dfs_traverse(node):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]
            depth = node.get("depth", 0)

            # Check for Boc protection reaction
            if (
                checker.check_reaction("Boc amine protection", rsmi)
                or checker.check_reaction("Boc amine protection explicit", rsmi)
                or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                or checker.check_reaction("Boc amine protection of primary amine", rsmi)
            ):

                print(f"Boc protection detected at depth: {depth}")
                # The product is now Boc-protected
                protected_molecules["boc"][product] = depth

            # Check for Boc deprotection reaction
            elif (
                checker.check_reaction("Boc amine deprotection", rsmi)
                or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
            ):

                print(f"Boc deprotection detected at depth: {depth}")
                # The product is now deprotected
                deprotected_molecules["boc"][product] = depth

            # Check for Cbz protection
            elif checker.check_reaction(
                "Carboxyl benzyl deprotection", rsmi
            ) or checker.check_reaction("Hydroxyl benzyl deprotection", rsmi):

                print(f"Cbz protection detected at depth: {depth}")
                protected_molecules["cbz"][product] = depth

            # Check for Cbz deprotection
            elif checker.check_fg("Carbamic ester", product) and any(
                checker.check_fg("Carbamic ester", r) for r in reactants
            ):
                # This is a simplified check for Cbz deprotection
                # Ideally, we would check for specific Cbz deprotection reactions
                print(f"Cbz deprotection detected at depth: {depth}")
                deprotected_molecules["cbz"][product] = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have a valid protection-deprotection sequence
    boc_sequence = False
    cbz_sequence = False

    # For Boc: Check if any molecule was protected and later deprotected
    for protected_mol, protection_depth in protected_molecules["boc"].items():
        for deprotected_mol, deprotection_depth in deprotected_molecules["boc"].items():
            # In retrosynthesis, protection happens at higher depth than deprotection
            if deprotection_depth < protection_depth:
                boc_sequence = True
                print(f"Complete Boc protection-deprotection sequence detected")
                break

    # For Cbz: Check if any molecule was protected and later deprotected
    for protected_mol, protection_depth in protected_molecules["cbz"].items():
        for deprotected_mol, deprotection_depth in deprotected_molecules["cbz"].items():
            # In retrosynthesis, protection happens at higher depth than deprotection
            if deprotection_depth < protection_depth:
                cbz_sequence = True
                print(f"Complete Cbz protection-deprotection sequence detected")
                break

    # If we have at least one deprotection event, that's also a sign of a protection-deprotection strategy
    boc_deprotection_events = len(deprotected_molecules["boc"])
    cbz_deprotection_events = len(deprotected_molecules["cbz"])

    # Determine if our strategy is present
    protection_strategy_present = (
        boc_sequence or cbz_sequence or boc_deprotection_events > 0 or cbz_deprotection_events > 0
    )

    print(f"Protection-deprotection strategy detected: {protection_strategy_present}")
    print(f"Boc protection events: {len(protected_molecules['boc'])}")
    print(f"Boc deprotection events: {boc_deprotection_events}")
    print(f"Cbz protection events: {len(protected_molecules['cbz'])}")
    print(f"Cbz deprotection events: {cbz_deprotection_events}")

    return protection_strategy_present
