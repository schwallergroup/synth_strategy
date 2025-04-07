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


def main(synthetic_route):
    """
    Detects if a reaction in the synthetic route is a Suzuki coupling.

    Args:
        synthetic_route: JSON object representing the synthetic route

    Returns:
        bool: True if any reaction in the route is a Suzuki coupling
    """

    def dfs_traverse(node):
        # If this is a reaction node, check if it's a Suzuki coupling
        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Checking reaction: {rsmi}")

            # Check if this is a Suzuki coupling using the checker function
            suzuki_types = [
                "Suzuki coupling with boronic acids",
                "Suzuki coupling with boronic acids OTf",
                "Suzuki coupling with sulfonic esters",
                "Suzuki coupling with boronic esters OTf",
                "Suzuki coupling with boronic esters",
                "{Suzuki}",
            ]

            for suzuki_type in suzuki_types:
                if checker.check_reaction(suzuki_type, rsmi):
                    print(f"Found Suzuki coupling: {suzuki_type}")
                    return True

            # If the reaction checker didn't identify it, try a more detailed check
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl halide in reactants
                has_aryl_halide = any(checker.check_fg("Aromatic halide", r) for r in reactants)

                # Check for boronic acid/ester in reactants
                has_boronic_acid = any(checker.check_fg("Boronic acid", r) for r in reactants)
                has_boronic_ester = any(checker.check_fg("Boronic ester", r) for r in reactants)

                if has_aryl_halide and (has_boronic_acid or has_boronic_ester):
                    print("Found Suzuki coupling based on functional group analysis")
                    return True
            except Exception as e:
                print(f"Error analyzing reaction components: {e}")

        # Recursively check children
        for child in node.get("children", []):
            if dfs_traverse(child):
                return True

        return False

    # Start traversal from the root
    return dfs_traverse(synthetic_route)
