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
    This function detects if the synthesis includes a late-stage imine formation.
    Late-stage is defined as occurring within the final 3 steps (depth <= 2).
    """
    imine_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal imine_formation_detected

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if (
            node["type"] == "reaction"
            and depth <= 2
            and "rsmi" in node.get("metadata", {})
        ):
            rsmi = node["metadata"]["rsmi"]
            print(f"Checking reaction at depth {depth}: {rsmi}")

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if product contains an imine
            product_has_imine = checker.check_fg(
                "Substituted imine", product_smiles
            ) or checker.check_fg("Unsubstituted imine", product_smiles)

            # Check if any reactant contains an imine
            reactants_have_imine = any(
                checker.check_fg("Substituted imine", r)
                or checker.check_fg("Unsubstituted imine", r)
                for r in reactants_smiles
            )

            print(
                f"Product has imine: {product_has_imine}, Reactants have imine: {reactants_have_imine}"
            )

            # Check for imine formation reactions
            imine_reaction = (
                checker.check_reaction(
                    "Addition of primary amines to aldehydes/thiocarbonyls", rsmi
                )
                or checker.check_reaction(
                    "Addition of primary amines to ketones/thiocarbonyls", rsmi
                )
                or checker.check_reaction(
                    "Addition of secondary amines to aldehydes/thiocarbonyls", rsmi
                )
                or checker.check_reaction(
                    "Addition of secondary amines to ketones/thiocarbonyls", rsmi
                )
                or checker.check_reaction("Ketone/aldehyde to hydrazone", rsmi)
                or checker.check_reaction("reductive amination", rsmi)
            )

            print(f"Imine formation reaction detected: {imine_reaction}")

            # Additional check for amine and carbonyl reactants
            has_amine = any(
                checker.check_fg("Primary amine", r)
                or checker.check_fg("Secondary amine", r)
                for r in reactants_smiles
            )

            has_carbonyl = any(
                checker.check_fg("Aldehyde", r)
                or checker.check_fg("Ketone", r)
                or checker.check_fg("Formaldehyde", r)
                for r in reactants_smiles
            )

            print(
                f"Reactants have amine: {has_amine}, Reactants have carbonyl: {has_carbonyl}"
            )

            # Detect imine formation: either by pattern change, specific reaction, or reactant combination
            if (
                (product_has_imine and not reactants_have_imine)
                or imine_reaction
                or (product_has_imine and has_amine and has_carbonyl)
            ):
                imine_formation_detected = True
                print(f"Late-stage imine formation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {imine_formation_detected}")
    return imine_formation_detected
