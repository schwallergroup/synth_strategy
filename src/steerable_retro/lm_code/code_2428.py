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
    Detects a synthetic strategy involving sequential functionalization of a bicyclic aromatic core.
    """
    # List of common bicyclic aromatic systems to check
    bicyclic_rings = [
        "naphthalene",
        "quinoline",
        "isoquinoline",
        "indole",
        "benzothiophene",
        "benzofuran",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "purine",
        "carbazole",
        "dibenzofuran",
        "dibenzothiophene",
        "acridine",
        "anthracene",
    ]

    # List of functional groups that could be added to the core
    functional_groups = [
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Nitrile",
        "Carboxylic acid",
        "Ester",
        "Amide",
        "Aldehyde",
        "Ketone",
        "Alcohol",
        "Ether",
        "Halide",
        "Nitro group",
        "Sulfonamide",
        "Boronic acid",
        "Boronic ester",
        "Alkyne",
        "Alkene",
    ]

    # Track modifications by depth
    modifications_by_depth = {}

    # Track core formation
    core_formation = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains a bicyclic aromatic core
                has_bicyclic_core = False
                core_type = None

                for ring in bicyclic_rings:
                    if checker.check_ring(ring, product):
                        has_bicyclic_core = True
                        core_type = ring
                        break

                if has_bicyclic_core:
                    # Check if the core exists in at least one reactant
                    core_in_reactants = False
                    for reactant in reactants:
                        if checker.check_ring(core_type, reactant):
                            core_in_reactants = True
                            break

                    # If core doesn't exist in reactants, this is a core-forming reaction
                    if not core_in_reactants:
                        print(f"Found core-forming reaction for {core_type} at depth {depth}")
                        core_formation[depth] = {"core": core_type, "reaction": rsmi}
                    else:
                        # Check for common functionalization reactions
                        functionalization_reactions = [
                            "Suzuki coupling",
                            "Buchwald-Hartwig",
                            "N-arylation",
                            "Heck",
                            "Sonogashira",
                            "Friedel-Crafts acylation",
                            "Friedel-Crafts alkylation",
                            "Aromatic chlorination",
                            "Aromatic bromination",
                            "Aromatic iodination",
                            "Aromatic fluorination",
                            "Aromatic nitration",
                            "Minisci (para)",
                            "Minisci (ortho)",
                            "Chan-Lam amine",
                            "Chan-Lam alcohol",
                            "Chan-Lam etherification",
                            "Negishi coupling",
                            "Stille reaction_aryl",
                            "Hiyama-Denmark Coupling",
                            "Ullmann-Goldberg Substitution amine",
                            "Ullmann-Goldberg Substitution thiol",
                            "Ullmann-Goldberg Substitution aryl alcohol",
                        ]

                        is_functionalization = False
                        reaction_type = None

                        # First check for known functionalization reactions
                        for rxn_type in functionalization_reactions:
                            if checker.check_reaction(rxn_type, rsmi):
                                is_functionalization = True
                                reaction_type = rxn_type
                                break

                        # If not a known reaction type, check for functional group changes
                        if not is_functionalization:
                            # Check if functional groups are added/modified on the product
                            product_fgs = []
                            for fg in functional_groups:
                                if checker.check_fg(fg, product):
                                    product_fgs.append(fg)

                            # Check if these functional groups are new or modified compared to reactants
                            for reactant in reactants:
                                if checker.check_ring(core_type, reactant):
                                    reactant_fgs = []
                                    for fg in functional_groups:
                                        if checker.check_fg(fg, reactant):
                                            reactant_fgs.append(fg)

                                    # If product has functional groups not in reactant, it's a modification
                                    new_fgs = [fg for fg in product_fgs if fg not in reactant_fgs]
                                    if new_fgs:
                                        is_functionalization = True
                                        reaction_type = f"Addition of {', '.join(new_fgs)}"
                                        break

                        if is_functionalization:
                            print(
                                f"Found bicyclic core ({core_type}) modification at depth {depth} via {reaction_type}"
                            )
                            modifications_by_depth[depth] = {
                                "core": core_type,
                                "reaction": reaction_type,
                            }

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have at least 2 modifications or 1 core formation + 1 modification
    total_modifications = len(modifications_by_depth)
    has_core_formation = len(core_formation) > 0

    print(f"Found {total_modifications} modifications and {len(core_formation)} core formations")

    if total_modifications < 2 and not (has_core_formation and total_modifications >= 1):
        print("Insufficient modifications found")
        return False

    # If we have core formation and modifications, check if they're sequential
    if has_core_formation and total_modifications >= 1:
        core_depths = sorted(core_formation.keys())
        mod_depths = sorted(modifications_by_depth.keys())

        # Check if core formation is followed by modification
        for core_depth in core_depths:
            for mod_depth in mod_depths:
                if mod_depth > core_depth and mod_depth - core_depth <= 3:
                    core_type_formation = core_formation[core_depth]["core"]
                    core_type_modification = modifications_by_depth[mod_depth]["core"]

                    if core_type_formation == core_type_modification:
                        print(
                            f"Found core formation at depth {core_depth} followed by modification at depth {mod_depth}"
                        )
                        return True

    # Check if all modifications are on the same core type
    if total_modifications >= 2:
        core_types = set(mod["core"] for mod in modifications_by_depth.values())
        if len(core_types) > 1:
            print("Modifications are on different core types")
            return False

        # Check if the modifications are sequential
        depths = sorted(modifications_by_depth.keys())

        # Check if at least two modifications are within 3 steps of each other
        for i in range(len(depths) - 1):
            if depths[i + 1] - depths[i] <= 3:  # Allow up to 3 steps between modifications
                print(f"Found sequential modifications at depths {depths[i]} and {depths[i+1]}")
                return True

    print("No sequential modifications found")
    return False
