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
    Detects if the route contains a Suzuki coupling in the second half of the synthesis.
    Looks for a reaction where an aryl halide and boronic acid form a biaryl.
    """
    suzuki_found = False
    max_depth = 0

    # First pass to find max_depth
    def find_max_depth(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)
        for child in node.get("children", []):
            find_max_depth(child, depth + 1)

    find_max_depth(route)
    print(f"Maximum synthesis depth: {max_depth}")

    # Second pass to find Suzuki couplings
    def dfs_traverse(node, depth=0):
        nonlocal suzuki_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Check if this is a Suzuki coupling directly
            is_suzuki = (
                checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with sulfonic esters", rsmi)
            )

            if is_suzuki:
                print(f"Found Suzuki coupling at depth {depth}")

                # Check if this is in the second half of the synthesis (late stage)
                if depth <= max_depth / 2:
                    print(f"This is a late-stage Suzuki coupling (depth {depth} <= {max_depth/2})")
                    suzuki_found = True
            else:
                # Fallback check for Suzuki coupling by looking at functional groups
                has_aryl_halide = False
                has_boronic = False

                for reactant in reactants:
                    if not reactant:
                        continue

                    if checker.check_fg("Aromatic halide", reactant):
                        has_aryl_halide = True
                        print(f"Found aryl halide in reactant: {reactant}")

                    if checker.check_fg("Boronic acid", reactant) or checker.check_fg(
                        "Boronic ester", reactant
                    ):
                        has_boronic = True
                        print(f"Found boronic acid/ester in reactant: {reactant}")

                if has_aryl_halide and has_boronic:
                    print(f"Found potential Suzuki coupling by functional groups at depth {depth}")

                    # Check if this is in the second half of the synthesis
                    if depth <= max_depth / 2:
                        print(
                            f"This is a late-stage potential Suzuki coupling (depth {depth} <= {max_depth/2})"
                        )
                        suzuki_found = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return suzuki_found
