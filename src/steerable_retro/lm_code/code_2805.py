#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    This function detects if the synthetic route involves protection of diols as cyclic acetals/ketals.
    """
    # Track if we found a valid protection reaction
    found_protection_reaction = False

    def dfs_traverse(node):
        nonlocal found_protection_reaction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction: {rsmi}")

            # Check for protection reactions (acetalization)
            if checker.check_reaction("Diol acetalization", rsmi) or checker.check_reaction(
                "Aldehyde or ketone acetalization", rsmi
            ):
                print(f"Found acetalization reaction: {rsmi}")

                # Check for cyclic acetal/ketal in product
                cyclic_acetal_rings = ["dioxolane", "dioxane", "dioxolene", "dioxepane", "trioxane"]
                has_cyclic_acetal = any(
                    checker.check_ring(ring, product) for ring in cyclic_acetal_rings
                )

                if has_cyclic_acetal:
                    print(f"Found cyclic acetal/ketal in product: {product}")
                    found_protection_reaction = True
                    return

            # Direct detection of cyclic acetal/ketal formation
            cyclic_acetal_rings = ["dioxolane", "dioxane", "dioxolene", "dioxepane", "trioxane"]
            has_cyclic_acetal_product = any(
                checker.check_ring(ring, product) for ring in cyclic_acetal_rings
            )
            has_cyclic_acetal_reactants = any(
                any(checker.check_ring(ring, reactant) for ring in cyclic_acetal_rings)
                for reactant in reactants
            )

            # If a cyclic acetal/ketal is formed (appears in product but not in reactants)
            if has_cyclic_acetal_product and not has_cyclic_acetal_reactants:
                print(f"Found cyclic acetal/ketal formation in product: {product}")

                # Check for diol in reactants
                has_diol = False
                has_carbonyl = False

                for reactant in reactants:
                    # Check for ethylene glycol or similar diols
                    if "OCO" in reactant or "OCCO" in reactant or "OCCCO" in reactant:
                        has_diol = True
                        print(f"Found potential diol pattern in reactant: {reactant}")

                    # Count alcohols in this specific reactant
                    alcohol_count = 0
                    if checker.check_fg("Primary alcohol", reactant):
                        alcohol_count += 1
                        print(f"Found primary alcohol in reactant: {reactant}")
                    if checker.check_fg("Secondary alcohol", reactant):
                        alcohol_count += 1
                        print(f"Found secondary alcohol in reactant: {reactant}")
                    if checker.check_fg("Tertiary alcohol", reactant):
                        alcohol_count += 1
                        print(f"Found tertiary alcohol in reactant: {reactant}")

                    # If this single reactant has multiple alcohols, it's a potential diol
                    if alcohol_count >= 2:
                        has_diol = True
                        print(f"Found potential diol in reactant: {reactant}")

                    # Check for carbonyl compounds
                    if (
                        checker.check_fg("Aldehyde", reactant)
                        or checker.check_fg("Ketone", reactant)
                        or checker.check_fg("Formaldehyde", reactant)
                    ):
                        has_carbonyl = True
                        print(f"Found carbonyl compound in reactant: {reactant}")

                # If we have either a diol or a carbonyl, and a cyclic acetal is formed, it's likely a protection
                if has_diol or has_carbonyl:
                    print("Found valid cyclic acetal/ketal formation")
                    found_protection_reaction = True
                    return

            # Check for deprotection reactions (hydrolysis)
            if (
                checker.check_reaction("Acetal hydrolysis to diol", rsmi)
                or checker.check_reaction("Acetal hydrolysis to aldehyde", rsmi)
                or checker.check_reaction("Ketal hydrolysis to ketone", rsmi)
            ):
                print(f"Found acetal/ketal hydrolysis reaction: {rsmi}")

                # Check for cyclic acetal/ketal in reactants
                cyclic_acetal_rings = ["dioxolane", "dioxane", "dioxolene", "dioxepane", "trioxane"]
                has_cyclic_acetal_reactant = any(
                    any(checker.check_ring(ring, reactant) for ring in cyclic_acetal_rings)
                    for reactant in reactants
                )

                if has_cyclic_acetal_reactant:
                    print("Found valid cyclic acetal/ketal deprotection")
                    found_protection_reaction = True
                    return

            # Special case: Check for reactions where a cyclic acetal/ketal appears or disappears
            # This catches cases where the reaction might not be explicitly labeled as acetalization/hydrolysis
            cyclic_acetal_rings = ["dioxolane", "dioxane", "dioxolene", "dioxepane", "trioxane"]
            has_cyclic_acetal_product = any(
                checker.check_ring(ring, product) for ring in cyclic_acetal_rings
            )
            has_cyclic_acetal_reactants = any(
                any(checker.check_ring(ring, reactant) for ring in cyclic_acetal_rings)
                for reactant in reactants
            )

            # Formation of cyclic acetal/ketal (protection)
            if has_cyclic_acetal_product and not has_cyclic_acetal_reactants:
                print(f"Detected formation of cyclic acetal/ketal: {product}")
                found_protection_reaction = True
                return

            # Removal of cyclic acetal/ketal (deprotection)
            if not has_cyclic_acetal_product and has_cyclic_acetal_reactants:
                print(f"Detected removal of cyclic acetal/ketal from reactants")
                found_protection_reaction = True
                return

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return found_protection_reaction
