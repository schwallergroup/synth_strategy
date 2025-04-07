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
    Detects a synthetic strategy involving late-stage epoxidation of a vinyl group.
    Looks for a vinyl group that gets converted to an oxirane (epoxide) ring in the final step.
    """
    vinyl_to_epoxide_found = False

    def dfs_traverse(node, depth=0):
        nonlocal vinyl_to_epoxide_found

        if node["type"] == "reaction" and depth <= 1:  # Final or penultimate reaction (late stage)
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                print(f"No reaction SMILES found in metadata at depth {depth}")
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if any reactant has a vinyl group but no epoxide
            has_vinyl_reactant = False
            vinyl_reactant = None

            for reactant in reactants_smiles:
                try:
                    if checker.check_fg("Vinyl", reactant) and not checker.check_ring(
                        "oxirane", reactant
                    ):
                        has_vinyl_reactant = True
                        vinyl_reactant = reactant
                        print(f"Found reactant with vinyl group: {reactant}")
                        break
                except Exception as e:
                    print(f"Error checking vinyl in reactant: {e}")
                    continue

            # Check if product has an epoxide ring
            has_epoxide_product = False

            try:
                if checker.check_ring("oxirane", product_smiles):
                    has_epoxide_product = True
                    print(f"Found product with epoxide ring: {product_smiles}")
            except Exception as e:
                print(f"Error checking epoxide in product: {e}")

            # Check if this is an epoxidation reaction
            is_epoxidation = False
            try:
                # Check for common epoxidation reaction types
                epoxidation_reactions = [
                    "Williamson Ether Synthesis (intra to epoxy)",
                    "Alkene to diol",
                    "anti-Markovnikov alkene hydration to alcohol",
                    "Markovnikov alkene hydration to alcohol",
                ]

                for rxn_type in epoxidation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_epoxidation = True
                        print(f"Identified as {rxn_type}")
                        break

                # If specific epoxidation reaction type not found, check for general pattern
                if not is_epoxidation and has_vinyl_reactant and has_epoxide_product:
                    # Verify the epoxide is newly formed (not in any reactant)
                    epoxide_in_reactants = False
                    for reactant in reactants_smiles:
                        if checker.check_ring("oxirane", reactant):
                            epoxide_in_reactants = True
                            print("Epoxide already present in reactant, not a new formation")
                            break

                    if not epoxide_in_reactants:
                        # Additional check: verify the reaction involves oxidation
                        # Common oxidants in epoxidation reactions
                        oxidants = ["O", "O2", "H2O2", "MCPBA", "mCPBA", "peracid", "peroxide"]
                        reagents = rsmi.split(">")[1]

                        if any(oxidant in reagents for oxidant in oxidants):
                            is_epoxidation = True
                            print("Identified as general epoxidation reaction (oxidation present)")
                        else:
                            # If no specific oxidant found, but we have vinyl → epoxide transformation
                            # It's likely still an epoxidation
                            is_epoxidation = True
                            print(
                                "Identified as general epoxidation reaction (vinyl → epoxide transformation)"
                            )
            except Exception as e:
                print(f"Error checking reaction type: {e}")

            # If all conditions are met, we have a vinyl to epoxide transformation
            if has_vinyl_reactant and has_epoxide_product and is_epoxidation:
                vinyl_to_epoxide_found = True
                print(f"Found epoxidation at depth {depth}: vinyl group converted to epoxide")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: vinyl_to_epoxide_found = {vinyl_to_epoxide_found}")

    return vinyl_to_epoxide_found
