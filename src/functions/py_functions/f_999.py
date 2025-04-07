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
    Detects if the synthesis involves formation of an epoxide from a hydroxyl group
    or other precursors, or epoxide ring opening reactions.
    """
    epoxide_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal epoxide_formation_found

        if node["type"] == "reaction" and not epoxide_formation_found:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if this is a known epoxide-forming reaction
                if checker.check_reaction(
                    "Williamson Ether Synthesis (intra to epoxy)", rsmi
                ):
                    print(
                        f"Epoxide formation detected via Williamson Ether Synthesis at depth {depth}"
                    )
                    epoxide_formation_found = True
                    return

                # If not a known reaction type, check manually
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Reactants: {reactants}")
                    print(f"Product: {product}")

                    # Check if any reactant has a hydroxyl group
                    reactant_has_hydroxyl = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        for r in reactants
                    )

                    # Check if any reactant has a halohydrin (another common epoxide precursor)
                    reactant_has_halohydrin = any(
                        (
                            checker.check_fg("Primary alcohol", r)
                            or checker.check_fg("Secondary alcohol", r)
                            or checker.check_fg("Tertiary alcohol", r)
                        )
                        and (
                            checker.check_fg("Primary halide", r)
                            or checker.check_fg("Secondary halide", r)
                            or checker.check_fg("Tertiary halide", r)
                        )
                        for r in reactants
                    )

                    # Check if any reactant has an epoxide ring
                    reactant_has_epoxide = any(
                        checker.check_ring("oxirane", r) for r in reactants
                    )

                    # Check if product has an epoxide ring
                    product_has_epoxide = checker.check_ring("oxirane", product)

                    # Check if product has a hydroxyl group (for epoxide opening)
                    product_has_hydroxyl = (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                    )

                    print(f"Reactant has hydroxyl: {reactant_has_hydroxyl}")
                    print(f"Reactant has halohydrin: {reactant_has_halohydrin}")
                    print(f"Reactant has epoxide: {reactant_has_epoxide}")
                    print(f"Product has epoxide: {product_has_epoxide}")
                    print(f"Product has hydroxyl: {product_has_hydroxyl}")

                    # Case 1: Epoxide formation - product has epoxide but reactants don't
                    if product_has_epoxide and not reactant_has_epoxide:
                        print(f"Epoxide formation detected at depth {depth}")
                        epoxide_formation_found = True
                        return

                    # Case 2: Epoxide opening - reactant has epoxide, product has hydroxyl
                    if (
                        reactant_has_epoxide
                        and product_has_hydroxyl
                        and not product_has_epoxide
                    ):
                        print(f"Epoxide ring opening detected at depth {depth}")
                        epoxide_formation_found = True
                        return

                    # Case 3: Hydroxyl to epoxide conversion
                    if (
                        reactant_has_hydroxyl or reactant_has_halohydrin
                    ) and product_has_epoxide:
                        print(
                            f"Hydroxyl to epoxide conversion detected at depth {depth}"
                        )
                        epoxide_formation_found = True
                        return

                except Exception as e:
                    print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: epoxide_formation_found = {epoxide_formation_found}")
    return epoxide_formation_found
