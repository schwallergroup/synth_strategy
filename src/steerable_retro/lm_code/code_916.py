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
    Detects Suzuki coupling reactions (aryl halide + boronic acid)
    """
    suzuki_coupling_found = False

    def dfs_traverse(node):
        nonlocal suzuki_coupling_found

        if node["type"] == "reaction" and not suzuki_coupling_found:
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction: {rsmi}")

                # Check for Suzuki coupling using the checker function
                suzuki_types = [
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic acids OTf",
                    "Suzuki coupling with sulfonic esters",
                    "Suzuki coupling with boronic esters OTf",
                    "Suzuki coupling with boronic esters",
                    "Suzuki",
                ]

                # First try to match known reaction types
                for suzuki_type in suzuki_types:
                    if checker.check_reaction(suzuki_type, rsmi):
                        suzuki_coupling_found = True
                        print(f"Found Suzuki coupling reaction: {suzuki_type}")
                        return

                # If no match found, check for characteristic functional groups
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactants contain boronic acid/ester and aryl halide
                has_boronic = False
                has_aryl_halide = False

                for reactant in reactants:
                    if checker.check_fg("Boronic acid", reactant) or checker.check_fg(
                        "Boronic ester", reactant
                    ):
                        has_boronic = True
                        print(f"Found boronic acid/ester in reactant: {reactant}")

                    if checker.check_fg("Aromatic halide", reactant):
                        has_aryl_halide = True
                        print(f"Found aromatic halide in reactant: {reactant}")

                # If both functional groups are present, it's likely a Suzuki coupling
                if has_boronic and has_aryl_halide:
                    suzuki_coupling_found = True
                    print("Found Suzuki coupling based on functional groups")
                else:
                    print("Not a Suzuki coupling reaction")

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        if not suzuki_coupling_found:
            for child in node.get("children", []):
                dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    print(f"Suzuki coupling found: {suzuki_coupling_found}")

    return suzuki_coupling_found
