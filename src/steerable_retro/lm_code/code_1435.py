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
    This function detects a combined strategy of sequential reductive aminations with
    Boc protection and thiazole modification.
    """
    # Initialize flags for each strategy
    reductive_amination_reactions = []
    boc_protection_reactions = []
    thiazole_reactions = []

    # Traverse the synthesis route to identify reactions
    def traverse_route(node, depth=0):
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]
                reactants = rsmi.split(">")[0].split(".")

                # Check for reductive amination
                if (
                    checker.check_reaction("Reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("Reductive amination with ketone", rsmi)
                    or checker.check_reaction("Reductive amination with alcohol", rsmi)
                    or checker.check_reaction("reductive amination", rsmi)
                ):
                    print(f"Found reductive amination at depth {depth}: {rsmi}")
                    reductive_amination_reactions.append((depth, rsmi))

                # Alternative check for reductive amination: look for amine formation
                elif any(
                    checker.check_fg("Primary amine", r)
                    or checker.check_fg("Secondary amine", r)
                    or checker.check_fg("Tertiary amine", r)
                    for r in reactants
                ) and (
                    checker.check_fg("Aldehyde", product) or checker.check_fg("Ketone", product)
                ):
                    print(
                        f"Found potential reductive amination (amine formation) at depth {depth}: {rsmi}"
                    )
                    reductive_amination_reactions.append((depth, rsmi))

                # Check for Boc protection
                if (
                    checker.check_reaction("Boc amine protection", rsmi)
                    or checker.check_reaction("Boc amine protection explicit", rsmi)
                    or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                    or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                    or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                    or checker.check_reaction("Boc amine protection of primary amine", rsmi)
                ):
                    print(f"Found Boc protection at depth {depth}: {rsmi}")
                    boc_protection_reactions.append((depth, rsmi))

                # Additional check for Boc group in product
                elif "OC(=O)C(C)(C)C" in product or "OC(O)=C(C)(C)C" in product:
                    print(f"Found Boc group in product at depth {depth}: {product}")
                    boc_protection_reactions.append((depth, rsmi))

                # Check for thiazole modification (focusing on thiazole and dechlorination)
                if checker.check_ring("thiazole", product):
                    print(f"Found thiazole in product at depth {depth}: {product}")

                    # Check if any reactant contains thiazole and chlorine
                    for reactant in reactants:
                        if checker.check_ring("thiazole", reactant) and checker.check_fg(
                            "Aromatic halide", reactant
                        ):
                            print(f"Found thiazole modification (dechlorination) at depth {depth}")
                            thiazole_reactions.append((depth, rsmi))
                            break
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Check for piperidine/piperazine in molecule nodes
        elif node["type"] == "mol" and depth == 0:  # Final product
            if checker.check_ring("piperidine", node["smiles"]) or checker.check_ring(
                "piperazine", node["smiles"]
            ):
                print(
                    f"Found piperidine/piperazine in final product, which often results from sequential reductive aminations"
                )

        # Recursively traverse children
        for child in node.get("children", []):
            traverse_route(child, depth + 1)

    # Start traversal from the root
    traverse_route(route)

    # Check for sequential reductive amination (at least 2 reductive amination reactions)
    has_sequential_reductive_amination = len(reductive_amination_reactions) >= 2

    # If no reductive aminations found, check if final product has piperidine/piperazine
    if not has_sequential_reductive_amination and route["type"] == "mol":
        if checker.check_ring("piperidine", route["smiles"]) or checker.check_ring(
            "piperazine", route["smiles"]
        ):
            print(
                f"Found piperidine/piperazine in final product, which often results from sequential reductive aminations"
            )
            has_sequential_reductive_amination = True

    # Check for Boc protection
    has_boc_protection = len(boc_protection_reactions) > 0

    # If no Boc protection reactions found, check if final product has Boc group
    if not has_boc_protection and route["type"] == "mol":
        if "OC(=O)C(C)(C)C" in route["smiles"] or "OC(O)=C(C)(C)C" in route["smiles"]:
            print(f"Found Boc group in final product")
            has_boc_protection = True

    # Check for thiazole modification
    has_thiazole_modification = len(thiazole_reactions) > 0

    print(f"Sequential reductive amination: {has_sequential_reductive_amination}")
    print(f"Boc protection: {has_boc_protection}")
    print(f"Thiazole modification: {has_thiazole_modification}")

    # Combined strategy requires at least reductive amination and one other feature
    if has_sequential_reductive_amination and (has_boc_protection or has_thiazole_modification):
        print("Detected combined reductive amination strategy with protection/modification")
        return True

    # Special case: If we have thiazole modification and the final product contains piperidine/piperazine with Boc
    if has_thiazole_modification and route["type"] == "mol":
        mol_smiles = route["smiles"]
        if (
            checker.check_ring("piperidine", mol_smiles)
            or checker.check_ring("piperazine", mol_smiles)
        ) and ("OC(=O)C(C)(C)C" in mol_smiles or "OC(O)=C(C)(C)C" in mol_smiles):
            print(
                "Detected combined strategy with thiazole modification and Boc-protected piperidine/piperazine"
            )
            return True

    return False
