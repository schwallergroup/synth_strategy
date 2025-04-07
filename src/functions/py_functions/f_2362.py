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
    This function detects if Boc protection is introduced early in the synthesis
    and maintained until later stages.
    """
    # Track Boc protection information
    boc_info = {
        "protected_molecules": [],  # List of molecules with Boc protection
        "protection_reactions": [],  # List of Boc protection reactions
        "deprotection_reactions": [],  # List of Boc deprotection reactions
        "depths": {},  # Dictionary mapping depths to Boc-protected molecules
    }

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check if molecule has Boc protection
            if checker.check_fg("Boc", mol_smiles):
                print(f"Found Boc-protected molecule at depth {depth}: {mol_smiles}")
                boc_info["protected_molecules"].append(mol_smiles)
                if depth not in boc_info["depths"]:
                    boc_info["depths"][depth] = []
                boc_info["depths"][depth].append(mol_smiles)

        elif (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for Boc protection reactions
            if (
                checker.check_reaction("Boc amine protection", rxn_smiles)
                or checker.check_reaction("Boc amine protection explicit", rxn_smiles)
                or checker.check_reaction(
                    "Boc amine protection with Boc anhydride", rxn_smiles
                )
                or checker.check_reaction(
                    "Boc amine protection (ethyl Boc)", rxn_smiles
                )
                or checker.check_reaction(
                    "Boc amine protection of secondary amine", rxn_smiles
                )
                or checker.check_reaction(
                    "Boc amine protection of primary amine", rxn_smiles
                )
            ):
                print(f"Found Boc protection reaction at depth {depth}: {rxn_smiles}")
                boc_info["protection_reactions"].append((depth, rxn_smiles))

            # Check for Boc deprotection reactions
            if (
                checker.check_reaction("Boc amine deprotection", rxn_smiles)
                or checker.check_reaction(
                    "Boc amine deprotection of guanidine", rxn_smiles
                )
                or checker.check_reaction(
                    "Boc amine deprotection to NH-NH2", rxn_smiles
                )
                or checker.check_reaction(
                    "Tert-butyl deprotection of amine", rxn_smiles
                )
            ):
                print(f"Found Boc deprotection reaction at depth {depth}: {rxn_smiles}")
                boc_info["deprotection_reactions"].append((depth, rxn_smiles))

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Traverse the route
    dfs_traverse(route)

    # Analyze the collected information
    depths = list(boc_info["depths"].keys())

    if not depths:
        print("No Boc-protected molecules found")
        return False

    # Check if Boc protection is introduced early (high depth)
    early_protection = max(depths) >= 3

    # Check if Boc protection is maintained until late stages (low depth)
    # Consider depth 2 or less as late stage
    late_stage_presence = min(depths) <= 2

    # Check if Boc is present in the final product (depth 0)
    boc_in_final_product = 0 in depths

    # Check if there are protection reactions
    has_protection_reactions = len(boc_info["protection_reactions"]) > 0

    # Check if Boc is present at multiple depths (maintained through synthesis)
    multiple_depths = len(depths) >= 2

    print(f"Boc protection depths: {depths}")
    print(f"Early protection: {early_protection}")
    print(f"Late stage presence: {late_stage_presence}")
    print(f"Boc in final product: {boc_in_final_product}")
    print(f"Has protection reactions: {has_protection_reactions}")
    print(f"Multiple depths: {multiple_depths}")

    # Return True if Boc is introduced early and maintained until late stages
    if early_protection and late_stage_presence and multiple_depths:
        print(
            "Boc protection strategy detected: introduced early and maintained until late stages"
        )
        return True

    return False
