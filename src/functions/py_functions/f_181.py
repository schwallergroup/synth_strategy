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
    This function detects if the synthesis route involves a Suzuki coupling reaction.
    """
    has_suzuki_coupling = False

    def dfs_traverse(node):
        nonlocal has_suzuki_coupling

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]

                # Use the checker function to directly check for Suzuki coupling reaction
                if (
                    checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction(
                        "Suzuki coupling with boronic acids OTf", rsmi
                    )
                    or checker.check_reaction(
                        "Suzuki coupling with boronic esters", rsmi
                    )
                    or checker.check_reaction(
                        "Suzuki coupling with boronic esters OTf", rsmi
                    )
                    or checker.check_reaction(
                        "Suzuki coupling with sulfonic esters", rsmi
                    )
                ):
                    has_suzuki_coupling = True
                    print(f"Detected Suzuki coupling: {rsmi}")
                else:
                    # Fallback to manual checking if the reaction checker doesn't identify it
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Check for boronic acid/ester in reactants
                    has_boronic = any(
                        checker.check_fg("Boronic acid", r)
                        or checker.check_fg("Boronic ester", r)
                        for r in reactants_smiles
                    )

                    # Check for aryl halide or triflate in reactants
                    has_electrophile = any(
                        checker.check_fg("Aromatic halide", r)
                        or checker.check_fg("Triflate", r)
                        for r in reactants_smiles
                    )

                    # If we have both key components, it might be a Suzuki coupling
                    if has_boronic and has_electrophile:
                        # Additional check: verify a new C-C bond is formed
                        # This is a simplified check - in a real implementation,
                        # we would need to track atom mappings to confirm
                        has_suzuki_coupling = True
                        print(
                            f"Detected potential Suzuki coupling (manual check): {rsmi}"
                        )
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Suzuki coupling strategy: {has_suzuki_coupling}")
    return has_suzuki_coupling
