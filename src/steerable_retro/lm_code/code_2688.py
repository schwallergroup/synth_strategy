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
    This function detects late-stage amide formation from carboxylic acid and amine.
    """
    late_amide_formation = False

    # List of amide formation reaction types to check
    amide_formation_reactions = [
        "Carboxylic acid with primary amine to amide",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Ester with primary amine to amide",
        "Ester with secondary amine to amide",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with secondary amine to amide",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal late_amide_formation

        if node["type"] == "reaction" and depth <= 2:  # Late stage (low depth)
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if this is an amide formation reaction
                    is_amide_formation = False
                    for reaction_type in amide_formation_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            is_amide_formation = True
                            print(f"Detected {reaction_type} reaction at depth {depth}")
                            break

                    if not is_amide_formation:
                        # Fallback: check for functional groups
                        has_carboxylic_acid = False
                        has_amine = False
                        has_amide_product = False

                        for reactant in reactants:
                            if checker.check_fg("Carboxylic acid", reactant):
                                has_carboxylic_acid = True
                                print(f"Found carboxylic acid in reactant: {reactant}")
                            if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                                "Secondary amine", reactant
                            ):
                                has_amine = True
                                print(f"Found amine in reactant: {reactant}")

                        if (
                            checker.check_fg("Primary amide", product)
                            or checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Tertiary amide", product)
                        ):
                            has_amide_product = True
                            print(f"Found amide in product: {product}")

                        if has_carboxylic_acid and has_amine and has_amide_product:
                            is_amide_formation = True
                            print(
                                f"Detected amide formation based on functional groups at depth {depth}"
                            )

                    if is_amide_formation:
                        print(f"Confirmed late-stage amide formation at depth {depth}")
                        late_amide_formation = True
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_amide_formation
