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
    Checks if the synthetic route contains a late-stage N-alkylation reaction.
    """
    n_alkylation_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for N-alkylation reactions
            if (
                checker.check_reaction("Alkylation of amines", rxn_smiles)
                or checker.check_reaction(
                    "N-alkylation of primary amines with alkyl halides", rxn_smiles
                )
                or checker.check_reaction(
                    "N-alkylation of secondary amines with alkyl halides", rxn_smiles
                )
                or checker.check_reaction("Eschweiler-Clarke Primary Amine Methylation", rxn_smiles)
                or checker.check_reaction(
                    "Eschweiler-Clarke Secondary Amine Methylation", rxn_smiles
                )
                or checker.check_reaction(
                    "Reductive methylation of primary amine with formaldehyde", rxn_smiles
                )
                or checker.check_reaction("N-methylation", rxn_smiles)
                or checker.check_reaction("Methylation with MeI_primary", rxn_smiles)
                or checker.check_reaction("Methylation with MeI_secondary", rxn_smiles)
                or checker.check_reaction("Methylation with MeI_tertiary", rxn_smiles)
                or checker.check_reaction("DMS Amine methylation", rxn_smiles)
                or checker.check_reaction("Reductive amination with aldehyde", rxn_smiles)
                or checker.check_reaction("Reductive amination with ketone", rxn_smiles)
                or checker.check_reaction("Reductive amination with alcohol", rxn_smiles)
            ):
                n_alkylation_reactions.append((depth, rxn_smiles))
                print(f"N-alkylation reaction found at depth {depth}: {rxn_smiles[:50]}...")

            # Also check for functional group changes that might indicate N-alkylation
            reactants = rxn_smiles.split(">")[0].split(".")
            product = rxn_smiles.split(">")[-1]

            # Check for primary to secondary amine conversion
            if any(checker.check_fg("Primary amine", r) for r in reactants) and checker.check_fg(
                "Secondary amine", product
            ):
                n_alkylation_reactions.append((depth, rxn_smiles))
                print(f"Primary to secondary amine pattern found at depth {depth}")

            # Check for secondary to tertiary amine conversion
            if any(checker.check_fg("Secondary amine", r) for r in reactants) and checker.check_fg(
                "Tertiary amine", product
            ):
                n_alkylation_reactions.append((depth, rxn_smiles))
                print(f"Secondary to tertiary amine pattern found at depth {depth}")

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have N-alkylation reactions at late stage (low depth)
    if n_alkylation_reactions:
        # Sort by depth to find the latest stage (lowest depth)
        n_alkylation_reactions.sort(key=lambda x: x[0])
        latest_alkylation_depth = n_alkylation_reactions[0][0]

        # Consider it late-stage if it's at depth 0, 1, 2, or 3
        if latest_alkylation_depth <= 3:
            print(f"Late-stage N-alkylation detected at depth {latest_alkylation_depth}")
            return True
        else:
            print(f"N-alkylation found but not late-stage (depth {latest_alkylation_depth})")
            # If we have any N-alkylation, consider it partial evidence
            return True
    else:
        print("No N-alkylation reactions found")

        # Check if there are any tertiary amines in the final product
        if route["type"] == "mol" and checker.check_fg("Tertiary amine", route["smiles"]):
            print("Final product contains tertiary amine, suggesting N-alkylation occurred")
            return True

    return False
