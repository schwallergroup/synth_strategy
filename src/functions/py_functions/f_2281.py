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
    This function detects if the synthesis employs a late-stage reductive amination
    to connect major fragments.
    """
    found_late_stage_reductive_amination = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_reductive_amination

        # For reaction nodes, check if it's a late-stage reductive amination
        if node["type"] == "reaction" and depth <= 2:  # Late stage (low depth)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Direct check for reductive amination reactions
                if (
                    checker.check_reaction("Reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("Reductive amination with ketone", rsmi)
                    or checker.check_reaction("Reductive amination with alcohol", rsmi)
                ):
                    print(
                        f"Found late-stage reductive amination at depth {depth}: {rsmi}"
                    )
                    found_late_stage_reductive_amination = True
                else:
                    # Fallback to functional group analysis
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for carbonyl and amine patterns in reactants
                    carbonyl_present = False
                    amine_present = False

                    for reactant in reactants:
                        if (
                            checker.check_fg("Aldehyde", reactant)
                            or checker.check_fg("Formaldehyde", reactant)
                            or checker.check_fg("Ketone", reactant)
                        ):
                            carbonyl_present = True
                            print(f"Found carbonyl group in reactant: {reactant}")

                        if checker.check_fg(
                            "Primary amine", reactant
                        ) or checker.check_fg("Secondary amine", reactant):
                            amine_present = True
                            print(f"Found amine group in reactant: {reactant}")

                    # Check if product has a new C-N bond (simplified check)
                    if carbonyl_present and amine_present:
                        # Check if the product has a new amine but no carbonyl
                        product_has_carbonyl = (
                            checker.check_fg("Aldehyde", product)
                            or checker.check_fg("Formaldehyde", product)
                            or checker.check_fg("Ketone", product)
                        )

                        product_has_amine = (
                            checker.check_fg("Primary amine", product)
                            or checker.check_fg("Secondary amine", product)
                            or checker.check_fg("Tertiary amine", product)
                        )

                        if product_has_amine and not product_has_carbonyl:
                            print(
                                f"Detected potential reductive amination at depth {depth}: {rsmi}"
                            )
                            found_late_stage_reductive_amination = True

        # Process children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Late-stage reductive amination detected: {found_late_stage_reductive_amination}"
    )
    return found_late_stage_reductive_amination
