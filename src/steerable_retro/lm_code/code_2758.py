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
    This function detects a synthetic strategy with a primarily linear sequence
    that includes one convergent Suzuki coupling step.
    """
    suzuki_coupling_found = False
    linear_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_coupling_found, linear_steps

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a Suzuki coupling using the checker function
            is_suzuki = (
                checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with sulfonic esters", rsmi)
                or checker.check_reaction("Suzuki", rsmi)
                or checker.check_reaction("{Suzuki}", rsmi)
            )

            # Manual check for Suzuki coupling patterns if checker fails
            if not is_suzuki:
                has_boronic = False
                has_aryl_halide = False

                for reactant in reactants:
                    # Check for boronic acid/ester
                    if "B(O)" in reactant or "BO" in reactant:
                        has_boronic = True
                        print(f"Found boronic acid/ester: {reactant}")

                    # Check for aryl halide
                    if any(x in reactant for x in ["Br", "I", "Cl"]) and any(
                        x in reactant for x in ["c", "C"]
                    ):
                        has_aryl_halide = True
                        print(f"Found aryl halide: {reactant}")

                # If both patterns are found, it's likely a Suzuki coupling
                if has_boronic and has_aryl_halide:
                    is_suzuki = True
                    print(f"Manually identified Suzuki coupling at depth {depth}")

            if is_suzuki:
                suzuki_coupling_found = True
                print(f"Suzuki coupling detected at depth {depth} with RSMI: {rsmi}")
            else:
                # Check if this is a linear step (only one non-starting material reactant)
                non_starting_material_count = 0
                for child in node.get("children", []):
                    if child["type"] == "mol" and not child.get("in_stock", False):
                        non_starting_material_count += 1

                if non_starting_material_count <= 1:
                    linear_steps += 1
                    print(f"Linear step detected at depth {depth} with RSMI: {rsmi}")
                else:
                    print(f"Non-linear step detected at depth {depth} with RSMI: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have a linear synthesis with one Suzuki coupling
    if suzuki_coupling_found and linear_steps >= 2:
        print(
            f"Strategy detected: Linear synthesis with one convergent Suzuki step ({linear_steps} linear steps)"
        )
        return True
    else:
        print(
            f"Strategy not detected: Suzuki coupling found: {suzuki_coupling_found}, Linear steps: {linear_steps}"
        )
        return False
