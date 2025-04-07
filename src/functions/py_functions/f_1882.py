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
    Detects if the route contains sequential carbonyl reductions (ester→alcohol, aldehyde→alcohol, ketone→alcohol)
    """
    carbonyl_reductions = []

    def dfs_traverse(node, depth=0):
        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Check for various carbonyl reductions
            is_reduction = False

            # Check for ester reduction
            if checker.check_reaction("Reduction of ester to primary alcohol", rsmi):
                for reactant in reactants:
                    if checker.check_fg("Ester", reactant) and checker.check_fg(
                        "Primary alcohol", product
                    ):
                        print(f"Found ester reduction at depth {depth}")
                        carbonyl_reductions.append((depth, "ester", rsmi))
                        is_reduction = True
                        break

            # Check for Grignard reactions with aldehydes (also a form of carbonyl reduction)
            elif checker.check_reaction("Grignard from aldehyde to alcohol", rsmi):
                for reactant in reactants:
                    if (
                        checker.check_fg("Aldehyde", reactant)
                        and checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                    ):
                        print(f"Found Grignard aldehyde reduction at depth {depth}")
                        carbonyl_reductions.append((depth, "aldehyde", rsmi))
                        is_reduction = True
                        break

            # Check for aldehyde/ketone reduction
            elif checker.check_reaction(
                "Reduction of aldehydes and ketones to alcohols", rsmi
            ):
                for reactant in reactants:
                    if checker.check_fg("Aldehyde", reactant) and checker.check_fg(
                        "Primary alcohol", product
                    ):
                        print(f"Found aldehyde reduction at depth {depth}")
                        carbonyl_reductions.append((depth, "aldehyde", rsmi))
                        is_reduction = True
                        break
                    elif checker.check_fg("Ketone", reactant) and checker.check_fg(
                        "Secondary alcohol", product
                    ):
                        print(f"Found ketone reduction at depth {depth}")
                        carbonyl_reductions.append((depth, "ketone", rsmi))
                        is_reduction = True
                        break

            # Check for carboxylic acid reduction
            elif checker.check_reaction(
                "Reduction of carboxylic acid to primary alcohol", rsmi
            ):
                for reactant in reactants:
                    if checker.check_fg(
                        "Carboxylic acid", reactant
                    ) and checker.check_fg("Primary alcohol", product):
                        print(f"Found carboxylic acid reduction at depth {depth}")
                        carbonyl_reductions.append((depth, "carboxylic acid", rsmi))
                        is_reduction = True
                        break

            # Additional reduction patterns
            elif checker.check_reaction("Reduction of nitrile to amine", rsmi):
                for reactant in reactants:
                    if checker.check_fg("Nitrile", reactant) and checker.check_fg(
                        "Primary amine", product
                    ):
                        print(f"Found nitrile reduction at depth {depth}")
                        carbonyl_reductions.append((depth, "nitrile", rsmi))
                        is_reduction = True
                        break

            # Check for amide reductions
            elif any(
                checker.check_reaction(rxn, rsmi)
                for rxn in [
                    "Reduction of primary amides to amines",
                    "Reduction of secondary amides to amines",
                    "Reduction of tertiary amides to amines",
                ]
            ):
                for reactant in reactants:
                    if (
                        checker.check_fg("Primary amide", reactant)
                        or checker.check_fg("Secondary amide", reactant)
                        or checker.check_fg("Tertiary amide", reactant)
                    ) and (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                    ):
                        print(f"Found amide reduction at depth {depth}")
                        carbonyl_reductions.append((depth, "amide", rsmi))
                        is_reduction = True
                        break

            # Check for functional group changes that indicate reductions
            if not is_reduction:
                # Check for ester to alcohol conversion (may not be caught by specific reaction check)
                for reactant in reactants:
                    if checker.check_fg("Ester", reactant) and checker.check_fg(
                        "Primary alcohol", product
                    ):
                        print(f"Found ester to alcohol conversion at depth {depth}")
                        carbonyl_reductions.append((depth, "ester", rsmi))
                        is_reduction = True
                        break

                # Check for aldehyde to alcohol conversion
                for reactant in reactants:
                    if checker.check_fg("Aldehyde", reactant) and (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                    ):
                        print(f"Found aldehyde to alcohol conversion at depth {depth}")
                        carbonyl_reductions.append((depth, "aldehyde", rsmi))
                        is_reduction = True
                        break

                # Check for ketone to alcohol conversion
                for reactant in reactants:
                    if checker.check_fg("Ketone", reactant) and checker.check_fg(
                        "Secondary alcohol", product
                    ):
                        print(f"Found ketone to alcohol conversion at depth {depth}")
                        carbonyl_reductions.append((depth, "ketone", rsmi))
                        is_reduction = True
                        break

                # Check for carboxylic acid to alcohol conversion
                for reactant in reactants:
                    if checker.check_fg(
                        "Carboxylic acid", reactant
                    ) and checker.check_fg("Primary alcohol", product):
                        print(
                            f"Found carboxylic acid to alcohol conversion at depth {depth}"
                        )
                        carbonyl_reductions.append((depth, "carboxylic acid", rsmi))
                        is_reduction = True
                        break

                # Check for nitrile to amine conversion
                for reactant in reactants:
                    if checker.check_fg("Nitrile", reactant) and checker.check_fg(
                        "Primary amine", product
                    ):
                        print(f"Found nitrile to amine conversion at depth {depth}")
                        carbonyl_reductions.append((depth, "nitrile", rsmi))
                        is_reduction = True
                        break

            # Check for oxidation of alcohol to aldehyde (reverse of reduction in retrosynthesis)
            if not is_reduction:
                for reactant in reactants:
                    if checker.check_fg(
                        "Primary alcohol", reactant
                    ) and checker.check_fg("Aldehyde", product):
                        print(
                            f"Found alcohol oxidation to aldehyde at depth {depth} (reverse reduction)"
                        )
                        carbonyl_reductions.append((depth, "alcohol_to_aldehyde", rsmi))
                        is_reduction = True
                        break
                    elif checker.check_fg(
                        "Secondary alcohol", reactant
                    ) and checker.check_fg("Ketone", product):
                        print(
                            f"Found alcohol oxidation to ketone at depth {depth} (reverse reduction)"
                        )
                        carbonyl_reductions.append((depth, "alcohol_to_ketone", rsmi))
                        is_reduction = True
                        break

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(
        f"Found {len(carbonyl_reductions)} carbonyl reductions: {carbonyl_reductions}"
    )

    # Sort reductions by depth
    carbonyl_reductions.sort(key=lambda x: x[0])

    # Check if we have at least 2 reductions and they are sequential (within 2 steps of each other)
    if len(carbonyl_reductions) >= 2:
        for i in range(len(carbonyl_reductions) - 1):
            if abs(carbonyl_reductions[i][0] - carbonyl_reductions[i + 1][0]) <= 2:
                print(
                    f"Found sequential carbonyl reductions: {carbonyl_reductions[i][1]} at depth {carbonyl_reductions[i][0]} and {carbonyl_reductions[i+1][1]} at depth {carbonyl_reductions[i+1][0]}"
                )
                return True

    return False
