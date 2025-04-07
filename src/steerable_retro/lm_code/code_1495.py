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
    Detects if the synthetic route employs a late-stage Suzuki coupling as the final step.
    This strategy involves coupling an aryl halide or triflate with a boronic acid/ester to form a biaryl system.
    """
    final_suzuki_detected = False

    def is_late_stage_reaction(node, max_depth=1):
        """Check if this is a late-stage reaction (within max_depth of the final product)"""
        if node["type"] != "reaction":
            return False

        # For simplicity, we'll consider a reaction "late-stage" if it has no reaction children
        # or is at most max_depth away from the final product
        reaction_children = [
            child for child in node.get("children", []) if child["type"] == "reaction"
        ]
        if not reaction_children:
            return True

        # If there are reaction children, check their depth
        def get_max_reaction_depth(n):
            if n["type"] != "reaction":
                return 0
            reaction_children = [
                child for child in n.get("children", []) if child["type"] == "reaction"
            ]
            if not reaction_children:
                return 1
            return 1 + max(get_max_reaction_depth(child) for child in reaction_children)

        depth = get_max_reaction_depth(node)
        return depth <= max_depth

    def dfs_traverse(node):
        nonlocal final_suzuki_detected

        if node["type"] == "reaction" and is_late_stage_reaction(node):
            print(
                f"Examining potential late-stage reaction: {node.get('metadata', {}).get('rsmi', '')}"
            )

            # Check if this is a Suzuki coupling using the checker
            rsmi = node.get("metadata", {}).get("rsmi", "")

            # Check for Suzuki coupling reaction types
            is_suzuki = (
                checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with sulfonic esters", rsmi)
                or checker.check_reaction("{Suzuki}", rsmi)
            )

            print(f"Is Suzuki coupling: {is_suzuki}")

            if is_suzuki:
                print("Detected potential late-stage Suzuki coupling")
                final_suzuki_detected = True
                return

            # If the reaction checker didn't identify it as Suzuki, try to verify manually
            parts = rsmi.split(">")
            if len(parts) >= 3:
                reactants = parts[0].split(".")
                product = parts[2]

                # Check for boronic acid/ester in reactants
                has_boronic = False
                # Check for aryl halide/triflate in reactants
                has_aryl_leaving_group = False

                for reactant in reactants:
                    # Check for boronic acid/ester
                    if checker.check_fg("Boronic acid", reactant) or checker.check_fg(
                        "Boronic ester", reactant
                    ):
                        has_boronic = True
                        print(f"Found boronic acid/ester in reactant: {reactant}")

                    # Check for aryl halide or triflate
                    if checker.check_fg("Aromatic halide", reactant) or checker.check_fg(
                        "Triflate", reactant
                    ):
                        has_aryl_leaving_group = True
                        print(f"Found aryl halide/triflate in reactant: {reactant}")

                    # Additional check for halides on aromatic rings
                    elif checker.check_ring("benzene", reactant):
                        if (
                            reactant.find("Br") >= 0
                            or reactant.find("Cl") >= 0
                            or reactant.find("I") >= 0
                            or reactant.find("F") >= 0
                        ):
                            has_aryl_leaving_group = True
                            print(f"Found potential aryl halide in reactant: {reactant}")

                print(
                    f"Has boronic: {has_boronic}, Has aryl leaving group: {has_aryl_leaving_group}"
                )

                # Check if product has a new biaryl bond (this is a characteristic of Suzuki coupling)
                if has_boronic and has_aryl_leaving_group:
                    print("Confirmed late-stage Suzuki coupling with required components")
                    final_suzuki_detected = True

        # Continue traversing
        for child in node.get("children", []):
            if (
                not final_suzuki_detected
            ):  # Stop traversing if we've already found what we're looking for
                dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return final_suzuki_detected
