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
    Detects the use of Boc protection in a synthetic route.

    A true Boc protection strategy involves:
    1. Protecting an amine with a Boc group
    2. Performing other chemistry while the amine is protected
    3. Deprotecting the Boc group later in the synthesis
    """
    # Track molecules with Boc protection
    boc_protected_molecules = set()
    # Track if we've found a complete protection-deprotection sequence
    found_complete_strategy = False
    # Track if we've found any Boc-related reactions
    found_any_boc_reaction = False

    def dfs_traverse(node, depth=0, path=None):
        nonlocal found_complete_strategy, found_any_boc_reaction

        if path is None:
            path = []

        current_path = path + [node]

        if node["type"] == "reaction":
            # Extract reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc protection reactions
            is_protection = any(
                checker.check_reaction(rxn, rsmi)
                for rxn in [
                    "Boc amine protection",
                    "Boc amine protection explicit",
                    "Boc amine protection with Boc anhydride",
                    "Boc amine protection (ethyl Boc)",
                    "Boc amine protection of secondary amine",
                    "Boc amine protection of primary amine",
                ]
            )

            # Check for Boc deprotection reactions
            is_deprotection = any(
                checker.check_reaction(rxn, rsmi)
                for rxn in [
                    "Boc amine deprotection",
                    "Boc amine deprotection of guanidine",
                    "Boc amine deprotection to NH-NH2",
                    "Tert-butyl deprotection of amine",
                ]
            )

            if is_protection:
                found_any_boc_reaction = True
                # Add the product to our set of Boc-protected molecules
                boc_protected_molecules.add(product)
                print(f"Found Boc protection reaction at depth {depth}")

            if is_deprotection:
                found_any_boc_reaction = True
                # Check if any of the reactants were previously protected
                for reactant in reactants:
                    if reactant in boc_protected_molecules:
                        found_complete_strategy = True
                        print(
                            f"Found complete Boc protection strategy: protected molecule was later deprotected"
                        )
                print(f"Found Boc deprotection reaction at depth {depth}")

            # Even if we don't have a direct protection-deprotection sequence,
            # check if the molecule contains a Boc group
            if not is_protection and not is_deprotection:
                for mol_smiles in reactants + [product]:
                    # Check if the molecule contains a Boc group
                    if checker.check_fg("Boc", mol_smiles):
                        found_any_boc_reaction = True
                        print(f"Found molecule with Boc group at depth {depth}")

        # For molecule nodes, check if they contain a Boc group
        elif node["type"] == "mol" and node["smiles"]:
            if checker.check_fg("Boc", node["smiles"]):
                found_any_boc_reaction = True
                print(f"Found molecule with Boc group at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found a complete protection-deprotection strategy
    # or if we found any Boc-related reactions or molecules
    if found_complete_strategy:
        print("Detected complete Boc protection-deprotection strategy")
        return True
    elif found_any_boc_reaction:
        print("Detected Boc protection/deprotection reactions or Boc groups")
        return True

    return False
