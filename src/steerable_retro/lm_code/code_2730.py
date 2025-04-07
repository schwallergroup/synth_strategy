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
    Detects if the route involves removal of one halogen and addition of another
    at different stages of the synthesis.
    """
    # Track halogen changes with depth information
    halogen_changes = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for direct halogen exchange reactions
                if checker.check_reaction("Finkelstein reaction", rsmi):
                    print(
                        f"Detected halogen exchange reaction: Finkelstein reaction at depth {depth}"
                    )
                    return True

                if checker.check_reaction("Aromatic substitution of bromine by chlorine", rsmi):
                    print(
                        f"Detected halogen exchange reaction: Aromatic substitution of bromine by chlorine at depth {depth}"
                    )
                    return True

                # Check for halogenation reactions
                halogenation_reactions = [
                    "Aromatic fluorination",
                    "Aromatic chlorination",
                    "Aromatic bromination",
                    "Aromatic iodination",
                    "Fluorination",
                    "Chlorination",
                    "Bromination",
                    "Iodination",
                ]

                for rxn_type in halogenation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected halogenation reaction: {rxn_type} at depth {depth}")
                        return True

                # Check for different halogen types in reactants
                reactant_halogens = {}
                for reactant in reactants:
                    if reactant.strip():
                        try:
                            # Check for various halogen functional groups
                            for halogen in ["I", "Br", "Cl", "F"]:
                                # Check primary halides
                                if checker.check_fg(f"Primary {halogen.lower()}ide", reactant):
                                    reactant_halogens[halogen] = (
                                        reactant_halogens.get(halogen, 0) + 1
                                    )

                                # Check secondary halides
                                if checker.check_fg(f"Secondary {halogen.lower()}ide", reactant):
                                    reactant_halogens[halogen] = (
                                        reactant_halogens.get(halogen, 0) + 1
                                    )

                                # Check tertiary halides
                                if checker.check_fg(f"Tertiary {halogen.lower()}ide", reactant):
                                    reactant_halogens[halogen] = (
                                        reactant_halogens.get(halogen, 0) + 1
                                    )

                                # Check aromatic halides
                                if checker.check_fg(f"Aromatic {halogen.lower()}ide", reactant):
                                    reactant_halogens[halogen] = (
                                        reactant_halogens.get(halogen, 0) + 1
                                    )

                                # Check alkenyl halides
                                if checker.check_fg(f"Alkenyl {halogen.lower()}ide", reactant):
                                    reactant_halogens[halogen] = (
                                        reactant_halogens.get(halogen, 0) + 1
                                    )

                                # Check haloalkynes
                                if halogen == "I" and checker.check_fg("Haloalkyne", reactant):
                                    reactant_halogens[halogen] = (
                                        reactant_halogens.get(halogen, 0) + 1
                                    )
                                elif halogen == "Br" and checker.check_fg("Haloalkyne", reactant):
                                    reactant_halogens[halogen] = (
                                        reactant_halogens.get(halogen, 0) + 1
                                    )
                                elif halogen == "Cl" and checker.check_fg("Haloalkyne", reactant):
                                    reactant_halogens[halogen] = (
                                        reactant_halogens.get(halogen, 0) + 1
                                    )
                                elif halogen == "F" and checker.check_fg("Haloalkyne", reactant):
                                    reactant_halogens[halogen] = (
                                        reactant_halogens.get(halogen, 0) + 1
                                    )
                        except Exception as e:
                            print(f"Error checking reactant halogens: {e}")
                            continue

                # Check for different halogen types in product
                product_halogens = {}
                try:
                    # Check for various halogen functional groups
                    for halogen in ["I", "Br", "Cl", "F"]:
                        # Check primary halides
                        if checker.check_fg(f"Primary {halogen.lower()}ide", product):
                            product_halogens[halogen] = product_halogens.get(halogen, 0) + 1

                        # Check secondary halides
                        if checker.check_fg(f"Secondary {halogen.lower()}ide", product):
                            product_halogens[halogen] = product_halogens.get(halogen, 0) + 1

                        # Check tertiary halides
                        if checker.check_fg(f"Tertiary {halogen.lower()}ide", product):
                            product_halogens[halogen] = product_halogens.get(halogen, 0) + 1

                        # Check aromatic halides
                        if checker.check_fg(f"Aromatic {halogen.lower()}ide", product):
                            product_halogens[halogen] = product_halogens.get(halogen, 0) + 1

                        # Check alkenyl halides
                        if checker.check_fg(f"Alkenyl {halogen.lower()}ide", product):
                            product_halogens[halogen] = product_halogens.get(halogen, 0) + 1

                        # Check haloalkynes
                        if halogen == "I" and checker.check_fg("Haloalkyne", product):
                            product_halogens[halogen] = product_halogens.get(halogen, 0) + 1
                        elif halogen == "Br" and checker.check_fg("Haloalkyne", product):
                            product_halogens[halogen] = product_halogens.get(halogen, 0) + 1
                        elif halogen == "Cl" and checker.check_fg("Haloalkyne", product):
                            product_halogens[halogen] = product_halogens.get(halogen, 0) + 1
                        elif halogen == "F" and checker.check_fg("Haloalkyne", product):
                            product_halogens[halogen] = product_halogens.get(halogen, 0) + 1
                except Exception as e:
                    print(f"Error checking product halogens: {e}")

                print(f"Reactant halogens at depth {depth}: {reactant_halogens}")
                print(f"Product halogens at depth {depth}: {product_halogens}")

                # Check for halogen exchange reactions
                for halogen in ["I", "Br", "Cl", "F"]:
                    reactant_count = reactant_halogens.get(halogen, 0)
                    product_count = product_halogens.get(halogen, 0)

                    if product_count < reactant_count:
                        # Halogen removal detected
                        halogen_changes.append(
                            {"type": "removal", "halogen": halogen, "depth": depth}
                        )
                        print(f"Detected {halogen} removal at depth {depth}")
                    elif product_count > reactant_count:
                        # Halogen addition detected
                        halogen_changes.append(
                            {"type": "addition", "halogen": halogen, "depth": depth}
                        )
                        print(f"Detected {halogen} addition at depth {depth}")

                # Check for direct exchange in a single reaction
                if len(reactant_halogens) > 0 and len(product_halogens) > 0:
                    for r_halogen, r_count in reactant_halogens.items():
                        for p_halogen, p_count in product_halogens.items():
                            if (
                                r_halogen != p_halogen
                                and r_count > product_halogens.get(r_halogen, 0)
                                and p_count > reactant_halogens.get(p_halogen, 0)
                            ):
                                print(
                                    f"Direct halogen exchange detected: {r_halogen} → {p_halogen} at depth {depth}"
                                )
                                return True

        # Continue traversing
        result = False
        for child in node.get("children", []):
            child_result = dfs_traverse(child, depth + 1)
            if child_result:
                result = True

        return result

    # Start traversal
    if dfs_traverse(route):
        return True

    # Check if we have both addition and removal of different halogens
    additions = [change for change in halogen_changes if change["type"] == "addition"]
    removals = [change for change in halogen_changes if change["type"] == "removal"]

    print(f"Total additions: {len(additions)}, Total removals: {len(removals)}")

    # If we have any halogen removals, consider it a halogen exchange strategy
    if len(removals) > 0:
        print(f"Halogen exchange strategy detected with {len(removals)} halogen removals")
        return True

    # Check for specific halogen exchange reactions that might have been missed
    def check_for_exchange_reactions(node):
        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]

            # Check for halogen exchange reactions
            exchange_reactions = [
                "Finkelstein reaction",
                "Aromatic substitution of bromine by chlorine",
                "Aromatic fluorination",
                "Aromatic chlorination",
                "Aromatic bromination",
                "Aromatic iodination",
                "Fluorination",
                "Chlorination",
                "Bromination",
                "Iodination",
            ]

            for rxn_type in exchange_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    print(f"Detected halogen exchange reaction: {rxn_type}")
                    return True

        # Continue traversing
        for child in node.get("children", []):
            if check_for_exchange_reactions(child):
                return True

        return False

    if check_for_exchange_reactions(route):
        return True

    return False
