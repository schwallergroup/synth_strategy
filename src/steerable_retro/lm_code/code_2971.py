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
    Detects if the route involves early-stage functional group manipulations
    (alcohol-aldehyde-ester interconversions) at depths > 1.
    """
    has_early_fg_manipulation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_early_fg_manipulation

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and depth > 1:  # Early stage (depth > 1)
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    print(f"No reaction SMILES found at depth {depth}")
                    return

                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")
                print(f"Reactants: {reactants}")
                print(f"Product: {product}")

                # Debug functional groups
                for i, r in enumerate(reactants):
                    print(
                        f"Reactant {i} FGs: Primary alcohol: {checker.check_fg('Primary alcohol', r)}, "
                        f"Secondary alcohol: {checker.check_fg('Secondary alcohol', r)}, "
                        f"Ester: {checker.check_fg('Ester', r)}, "
                        f"Aldehyde: {checker.check_fg('Aldehyde', r)}, "
                        f"Carboxylic acid: {checker.check_fg('Carboxylic acid', r)}"
                    )
                print(
                    f"Product FGs: Primary alcohol: {checker.check_fg('Primary alcohol', product)}, "
                    f"Secondary alcohol: {checker.check_fg('Secondary alcohol', product)}, "
                    f"Ester: {checker.check_fg('Ester', product)}, "
                    f"Aldehyde: {checker.check_fg('Aldehyde', product)}, "
                    f"Carboxylic acid: {checker.check_fg('Carboxylic acid', product)}"
                )

                # Check for alcohol-aldehyde-ester-carboxylic acid interconversions

                # Alcohol to aldehyde oxidation
                if any(
                    checker.check_fg("Primary alcohol", r) for r in reactants
                ) and checker.check_fg("Aldehyde", product):
                    print(f"Found alcohol to aldehyde oxidation at depth {depth}")
                    has_early_fg_manipulation = True

                # Aldehyde to alcohol reduction
                elif any(checker.check_fg("Aldehyde", r) for r in reactants) and checker.check_fg(
                    "Primary alcohol", product
                ):
                    print(f"Found aldehyde to alcohol reduction at depth {depth}")
                    has_early_fg_manipulation = True

                # Alcohol to ester transformation
                elif any(
                    checker.check_fg("Primary alcohol", r) for r in reactants
                ) and checker.check_fg("Ester", product):
                    print(f"Found alcohol to ester transformation at depth {depth}")
                    has_early_fg_manipulation = True

                # Ester to alcohol reduction
                elif any(checker.check_fg("Ester", r) for r in reactants) and checker.check_fg(
                    "Primary alcohol", product
                ):
                    print(f"Found ester to alcohol reduction at depth {depth}")
                    has_early_fg_manipulation = True

                # Carboxylic acid to ester
                elif any(
                    checker.check_fg("Carboxylic acid", r) for r in reactants
                ) and checker.check_fg("Ester", product):
                    print(f"Found carboxylic acid to ester transformation at depth {depth}")
                    has_early_fg_manipulation = True

                # Ester to carboxylic acid hydrolysis
                elif any(checker.check_fg("Ester", r) for r in reactants) and checker.check_fg(
                    "Carboxylic acid", product
                ):
                    print(f"Found ester to carboxylic acid hydrolysis at depth {depth}")
                    has_early_fg_manipulation = True

                # Aldehyde to carboxylic acid oxidation
                elif any(checker.check_fg("Aldehyde", r) for r in reactants) and checker.check_fg(
                    "Carboxylic acid", product
                ):
                    print(f"Found aldehyde to carboxylic acid oxidation at depth {depth}")
                    has_early_fg_manipulation = True

                # Secondary alcohol to ketone oxidation
                elif any(
                    checker.check_fg("Secondary alcohol", r) for r in reactants
                ) and checker.check_fg("Ketone", product):
                    print(f"Found secondary alcohol to ketone oxidation at depth {depth}")
                    has_early_fg_manipulation = True

                # Ketone to secondary alcohol reduction
                elif any(checker.check_fg("Ketone", r) for r in reactants) and checker.check_fg(
                    "Secondary alcohol", product
                ):
                    print(f"Found ketone to secondary alcohol reduction at depth {depth}")
                    has_early_fg_manipulation = True

                # Alcohol to carboxylic acid oxidation
                elif any(
                    checker.check_fg("Primary alcohol", r) for r in reactants
                ) and checker.check_fg("Carboxylic acid", product):
                    print(f"Found alcohol to carboxylic acid oxidation at depth {depth}")
                    has_early_fg_manipulation = True

                # Carboxylic acid to alcohol reduction
                elif any(
                    checker.check_fg("Carboxylic acid", r) for r in reactants
                ) and checker.check_fg("Primary alcohol", product):
                    print(f"Found carboxylic acid to alcohol reduction at depth {depth}")
                    has_early_fg_manipulation = True

                # Ketone to carboxylic acid oxidation
                elif any(checker.check_fg("Ketone", r) for r in reactants) and checker.check_fg(
                    "Carboxylic acid", product
                ):
                    print(f"Found ketone to carboxylic acid oxidation at depth {depth}")
                    has_early_fg_manipulation = True

                # Nitrile to carboxylic acid hydrolysis
                elif any(checker.check_fg("Nitrile", r) for r in reactants) and checker.check_fg(
                    "Carboxylic acid", product
                ):
                    print(f"Found nitrile to carboxylic acid hydrolysis at depth {depth}")
                    has_early_fg_manipulation = True

                # Nitrile to amide hydrolysis
                elif any(checker.check_fg("Nitrile", r) for r in reactants) and (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                ):
                    print(f"Found nitrile to amide hydrolysis at depth {depth}")
                    has_early_fg_manipulation = True

                # Amide to carboxylic acid hydrolysis
                elif (
                    any(checker.check_fg("Primary amide", r) for r in reactants)
                    or any(checker.check_fg("Secondary amide", r) for r in reactants)
                    or any(checker.check_fg("Tertiary amide", r) for r in reactants)
                ) and checker.check_fg("Carboxylic acid", product):
                    print(f"Found amide to carboxylic acid hydrolysis at depth {depth}")
                    has_early_fg_manipulation = True

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {has_early_fg_manipulation}")
    return has_early_fg_manipulation
