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
    This function detects late-stage amide formation (depth ≤ 2).
    It looks for a reaction where a carboxylic acid and an amine form an amide bond.
    """
    late_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_amide_formation

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a late-stage reaction (depth ≤ 2)
            if depth <= 2:
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check for amide formation reactions directly
                amide_reactions = [
                    "Carboxylic acid with primary amine to amide",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                ]

                for reaction_type in amide_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found late-stage amide formation: {reaction_type}")
                        late_amide_formation = True
                        return

                # If no specific reaction type matched, check for functional group changes
                has_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants)
                has_acyl_halide = any(checker.check_fg("Acyl halide", r) for r in reactants)
                has_ester = any(checker.check_fg("Ester", r) for r in reactants)

                has_primary_amine = any(checker.check_fg("Primary amine", r) for r in reactants)
                has_secondary_amine = any(checker.check_fg("Secondary amine", r) for r in reactants)

                has_primary_amide = checker.check_fg("Primary amide", product)
                has_secondary_amide = checker.check_fg("Secondary amide", product)
                has_tertiary_amide = checker.check_fg("Tertiary amide", product)

                # Check if we have the right reactants and product for amide formation
                if (
                    (has_acid or has_acyl_halide or has_ester)
                    and (has_primary_amine or has_secondary_amine)
                    and (has_primary_amide or has_secondary_amide or has_tertiary_amide)
                ):
                    print(f"Found late-stage amide formation through functional group analysis")
                    print(f"Acid: {has_acid}, Acyl halide: {has_acyl_halide}, Ester: {has_ester}")
                    print(
                        f"Primary amine: {has_primary_amine}, Secondary amine: {has_secondary_amine}"
                    )
                    print(
                        f"Primary amide: {has_primary_amide}, Secondary amide: {has_secondary_amide}, Tertiary amide: {has_tertiary_amide}"
                    )
                    late_amide_formation = True

        for child in node.get("children", []):
            if (
                not late_amide_formation
            ):  # Stop traversal if we already found what we're looking for
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_amide_formation
