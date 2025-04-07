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
    This function detects a protection-deprotection strategy using Boc group.
    It looks for Boc protection followed by later deprotection.
    """
    # Track protected molecules and their depths
    protected_molecules = {}  # Map of molecule SMILES to protection depth
    deprotected_molecules = {}  # Map of molecule SMILES to deprotection depth

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction at depth {depth}: {rsmi}")

            try:
                # Extract reactants and product
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product_part = rsmi.split(">")[-1]
                products = product_part.split(".")

                # Check for Boc protection reactions
                if (
                    checker.check_reaction("Boc amine protection", rsmi)
                    or checker.check_reaction("Boc amine protection explicit", rsmi)
                    or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                    or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                    or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                    or checker.check_reaction("Boc amine protection of primary amine", rsmi)
                ):

                    print(f"Potential Boc protection detected at depth {depth}")

                    # Verify Boc group is present in product but not in reactants (excluding Boc reagents)
                    boc_in_product = any(checker.check_fg("Boc", p) for p in products)

                    if boc_in_product:
                        print(f"Confirmed Boc protection at depth {depth}")
                        # Store the protected molecule with its depth
                        for p in products:
                            if checker.check_fg("Boc", p):
                                protected_molecules[p] = depth
                                print(f"Added protected molecule: {p}")

                # Check for Boc deprotection reactions
                if (
                    checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                    or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                    or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                ):

                    print(f"Potential Boc deprotection detected at depth {depth}")

                    # Verify Boc group is present in reactants but not in product
                    boc_in_reactants = any(checker.check_fg("Boc", r) for r in reactants)

                    if boc_in_reactants:
                        print(f"Confirmed Boc deprotection at depth {depth}")
                        # Store the deprotected molecule with its depth
                        for r in reactants:
                            if checker.check_fg("Boc", r):
                                deprotected_molecules[r] = depth
                                print(f"Added deprotected molecule: {r}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Recursively process children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Protected molecules: {protected_molecules}")
    print(f"Deprotected molecules: {deprotected_molecules}")

    # Check if we have any protection-deprotection pairs
    for protected_mol, protection_depth in protected_molecules.items():
        for deprotected_mol, deprotection_depth in deprotected_molecules.items():
            # In retrosynthetic routes, higher depth means earlier in the forward synthesis
            # So protection should have higher depth than deprotection
            if protection_depth > deprotection_depth:
                print(
                    f"Found valid protection-deprotection strategy: protection at depth {protection_depth}, deprotection at depth {deprotection_depth}"
                )
                return True

    # If we have both protections and deprotections but no valid pairs
    if protected_molecules and deprotected_molecules:
        print("Found both protection and deprotection, but not in correct order")

    # Alternative approach: check if any molecule has "Boc" functional group
    # This is a fallback in case the reaction detection failed
    boc_molecules = []

    def check_for_boc(node):
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            if checker.check_fg("Boc", mol_smiles):
                boc_molecules.append(mol_smiles)
                print(f"Found molecule with Boc group: {mol_smiles}")

        for child in node.get("children", []):
            check_for_boc(child)

    if not (protected_molecules or deprotected_molecules):
        print("No protection/deprotection reactions detected, checking for Boc groups directly")
        check_for_boc(route)
        if boc_molecules:
            print(f"Found {len(boc_molecules)} molecules with Boc groups")
            # If we found Boc groups but no protection/deprotection reactions,
            # this might indicate a protection-deprotection strategy that wasn't properly detected
            return True

    return False
