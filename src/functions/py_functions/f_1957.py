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
    This function detects if the synthesis includes a late-stage alcohol to halide conversion.
    """
    found_alcohol_to_halide = False

    # List of alcohol to halide conversion reaction types
    alcohol_to_halide_reactions = [
        "Alcohol to chloride_sulfonyl chloride",
        "Alcohol to chloride_SOCl2",
        "Alcohol to chloride_CHCl3",
        "Alcohol to chloride_CH2Cl2",
        "Alcohol to chloride_PCl5_ortho",
        "Alcohol to chloride_POCl3_ortho",
        "Alcohol to chloride_POCl3_para",
        "Alcohol to chloride_POCl3",
        "Alcohol to chloride_HCl",
        "Alcohol to chloride_Salt",
        "Alcohol to chloride_Other",
        "Appel reaction",  # For bromides and iodides
        "Primary alkyl halide to alcohol",  # Reverse reaction
        "Secondary alkyl halide to alcohol",  # Reverse reaction
    ]

    # List of alcohol functional groups
    alcohol_fgs = [
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Aromatic alcohol",
    ]

    # List of halide functional groups
    halide_fgs = [
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Aromatic halide",
        "Alkenyl halide",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_alcohol_to_halide

        if found_alcohol_to_halide:
            return

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a late-stage reaction (depth ≤ 2)
            if depth <= 2:
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Method 1: Check if this is a known alcohol-to-halide reaction type
                for reaction_type in alcohol_to_halide_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found {reaction_type} reaction at depth {depth}")
                        found_alcohol_to_halide = True
                        return

                # Method 2: Check for alcohol in reactants and halide in product
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for alcohols in reactants
                has_alcohol = False
                for reactant in reactants:
                    for alcohol_fg in alcohol_fgs:
                        if checker.check_fg(alcohol_fg, reactant):
                            print(f"{alcohol_fg} found in reactant: {reactant}")
                            has_alcohol = True
                            break
                    if has_alcohol:
                        break

                # Check for halides in product
                has_halide = False
                for halide_fg in halide_fgs:
                    if checker.check_fg(halide_fg, product):
                        print(f"{halide_fg} found in product: {product}")
                        has_halide = True
                        break

                # If both alcohol in reactants and halide in product, verify it's a conversion
                if has_alcohol and has_halide:
                    # Additional check: Make sure the halide wasn't already present in all reactants
                    all_reactants_have_halide = True
                    for reactant in reactants:
                        has_this_reactant_halide = False
                        for halide_fg in halide_fgs:
                            if checker.check_fg(halide_fg, reactant):
                                has_this_reactant_halide = True
                                break
                        if not has_this_reactant_halide:
                            all_reactants_have_halide = False
                            break

                    if not all_reactants_have_halide:
                        print(f"Alcohol to halide conversion detected at depth {depth}")
                        found_alcohol_to_halide = True
                        return

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_alcohol_to_halide
