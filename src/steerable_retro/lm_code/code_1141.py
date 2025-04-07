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
    Detects if the synthesis route involves a Suzuki coupling reaction.
    """
    suzuki_coupling_found = False

    def dfs_traverse(node):
        nonlocal suzuki_coupling_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a Suzuki coupling reaction using the checker function
            if (
                checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with sulfonic esters", rsmi)
                or checker.check_reaction("Suzuki", rsmi)
            ):

                print(f"Found Suzuki coupling reaction: {rsmi}")
                suzuki_coupling_found = True

            # Fallback check using functional groups if specific reaction check fails
            if not suzuki_coupling_found:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                has_boronic_acid = any(checker.check_fg("Boronic acid", r) for r in reactants)
                has_boronic_ester = any(checker.check_fg("Boronic ester", r) for r in reactants)
                has_aryl_halide = any(checker.check_fg("Aromatic halide", r) for r in reactants)

                if (has_boronic_acid or has_boronic_ester) and has_aryl_halide:
                    # Additional check for biaryl product formation would be ideal here
                    print(f"Found potential Suzuki coupling based on functional groups: {rsmi}")
                    suzuki_coupling_found = True

        # Recursively process children
        for child in node.get("children", []):
            if not suzuki_coupling_found:  # Optimization: stop traversal once found
                dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_coupling_found
