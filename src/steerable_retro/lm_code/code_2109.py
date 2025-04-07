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
    Detects a synthesis route that includes Boc protection of a primary amine.
    """
    boc_protection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_protection_found

        # Process reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Depth {depth}, Examining reaction: {rsmi}")

            # Check for Boc protection reactions (forward direction)
            boc_protection_reactions = [
                "Boc amine protection",
                "Boc amine protection explicit",
                "Boc amine protection with Boc anhydride",
                "Boc amine protection (ethyl Boc)",
                "Boc amine protection of primary amine",
                "Boc amine protection of secondary amine",
            ]

            # Check for Boc deprotection reactions (since we're traversing retrosynthetically)
            boc_deprotection_reactions = [
                "Boc amine deprotection",
                "Boc amine deprotection of guanidine",
                "Boc amine deprotection to NH-NH2",
                "Tert-butyl deprotection of amine",
            ]

            # Forward direction: Primary amine → Boc-protected amine
            if any(checker.check_reaction(rxn_type, rsmi) for rxn_type in boc_protection_reactions):
                # Check if any reactant has a primary or secondary amine
                amine_in_reactants = any(
                    checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                    for r in reactants
                )

                # Verify that the product has a Boc-protected amine
                has_boc_group = checker.check_fg("Carbamic ester", product)

                if amine_in_reactants and has_boc_group:
                    print(f"Found Boc protection reaction: {rsmi}")
                    boc_protection_found = True

            # Retrosynthetic direction: Boc-protected amine → Primary amine
            # This means we're looking at a deprotection reaction in the forward direction
            elif any(
                checker.check_reaction(rxn_type, rsmi) for rxn_type in boc_deprotection_reactions
            ):
                # Check if any reactant has a Boc-protected amine
                boc_in_reactants = any(checker.check_fg("Carbamic ester", r) for r in reactants)

                # Verify that the product has a primary amine
                has_primary_amine = checker.check_fg("Primary amine", product) or checker.check_fg(
                    "Secondary amine", product
                )

                if boc_in_reactants and has_primary_amine:
                    print(
                        f"Found Boc deprotection reaction (retrosynthetically implies protection): {rsmi}"
                    )
                    boc_protection_found = True

            # Special case: Check for Boc anhydride reaction
            elif any("CC(C)(C)OC(=O)OC(=O)OC(C)(C)C" in r for r in reactants):
                # Check if any reactant has a primary or secondary amine
                amine_in_reactants = any(
                    checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                    for r in reactants
                )

                # Verify that the product has a Boc-protected amine
                has_boc_group = checker.check_fg("Carbamic ester", product)

                if amine_in_reactants and has_boc_group:
                    print(f"Found Boc protection with Boc anhydride: {rsmi}")
                    boc_protection_found = True

            # Additional check for any reaction that produces a Boc-protected amine
            else:
                # Check if any reactant contains Boc anhydride or similar structure
                boc_reagent_in_reactants = any(
                    "OC(=O)OC(=O)O" in r or "C(C)(C)OC(=O)O" in r for r in reactants
                )

                # Check if any reactant has a primary or secondary amine
                amine_in_reactants = any(
                    checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                    for r in reactants
                )

                # Verify that the product has a Boc-protected amine
                has_boc_group = checker.check_fg("Carbamic ester", product)

                if boc_reagent_in_reactants and amine_in_reactants and has_boc_group:
                    print(f"Found Boc protection reaction (generic): {rsmi}")
                    boc_protection_found = True

        # Process molecule nodes (for debugging)
        elif node["type"] == "mol":
            smiles = node["smiles"]
            print(f"Depth {depth}, Examining molecule: {smiles}")
            if checker.check_fg("Primary amine", smiles):
                print(f"  - Contains primary amine")
            if checker.check_fg("Secondary amine", smiles):
                print(f"  - Contains secondary amine")
            if checker.check_fg("Carbamic ester", smiles):
                print(f"  - Contains Boc group (carbamic ester)")

        # Process children (retrosynthetic traversal)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Boc protection of amine found: {boc_protection_found}")
    return boc_protection_found
