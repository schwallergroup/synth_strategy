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
    Detects if the final step (depth 0 or 1) involves epoxide opening with an amine.
    """
    found_epoxide_opening = False

    def dfs_traverse(node, depth=0):
        nonlocal found_epoxide_opening

        # Update node depth
        node["depth"] = depth

        if node["type"] == "reaction" and depth <= 1:  # Check depth 0 and 1
            try:
                if "metadata" in node and "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]
                    print(f"Checking reaction at depth {depth}: {rsmi}")

                    # Check if this is an epoxide opening reaction with amine
                    if checker.check_reaction("Ring opening of epoxide with amine", rsmi):
                        print(f"Found late-stage epoxide opening with amine: {rsmi}")
                        found_epoxide_opening = True
                        return

                    # Fallback method: check reactants and products manually
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for epoxide in reactants
                    has_epoxide = any(checker.check_ring("oxirane", r) for r in reactants)
                    if has_epoxide:
                        print(f"Found epoxide in reactants")

                    # Check for amine in reactants
                    has_primary_amine = any(checker.check_fg("Primary amine", r) for r in reactants)
                    has_secondary_amine = any(
                        checker.check_fg("Secondary amine", r) for r in reactants
                    )
                    has_amine = has_primary_amine or has_secondary_amine

                    if has_amine:
                        print(
                            f"Found amine in reactants: Primary={has_primary_amine}, Secondary={has_secondary_amine}"
                        )

                    # Check if product has hydroxyl group (indicating epoxide opening)
                    has_hydroxyl = (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                        or checker.check_fg("Aromatic alcohol", product)
                    )

                    if has_hydroxyl:
                        print(f"Found hydroxyl group in product")

                    # Check if product has nitrogen-containing groups
                    has_n_product = (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                    )

                    if has_n_product:
                        print(f"Found nitrogen group in product")

                    # Determine if this is an epoxide opening reaction
                    if has_epoxide and has_amine and has_hydroxyl and has_n_product:
                        print(f"Found late-stage epoxide opening with amine (manual check): {rsmi}")
                        found_epoxide_opening = True

                    # Additional check for epoxide disappearing and amine reacting
                    if has_epoxide and has_amine:
                        # Check if epoxide is gone in product
                        if not checker.check_ring("oxirane", product):
                            print(f"Epoxide ring is opened in the reaction")
                            # If we have either hydroxyl or nitrogen in product, likely an epoxide opening
                            if has_hydroxyl or has_n_product:
                                print(
                                    f"Found late-stage epoxide opening with amine (ring disappearance): {rsmi}"
                                )
                                found_epoxide_opening = True
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_epoxide_opening
