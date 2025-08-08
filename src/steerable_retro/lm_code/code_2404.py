#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    Detects a synthetic strategy that uses nitrile functional groups as key intermediates.
    """
    has_nitrile_formation = False
    has_nitrile_transformation = False
    molecules_with_nitrile = []

    def dfs_traverse(node, depth=0):
        nonlocal has_nitrile_formation, has_nitrile_transformation

        # Check for nitrile in molecule nodes
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            if checker.check_fg("Nitrile", mol_smiles):
                molecules_with_nitrile.append((mol_smiles, depth))
                print(f"Found molecule with nitrile at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for nitrile presence in reactants and product
                reactants_have_nitrile = any(
                    checker.check_fg("Nitrile", r) for r in reactants_smiles if r
                )
                product_has_nitrile = (
                    checker.check_fg("Nitrile", product_smiles) if product_smiles else False
                )

                # In retrosynthetic traversal:
                # - Nitrile formation in forward direction appears as nitrile in product but not in reactants
                # - Nitrile transformation in forward direction appears as nitrile in reactants but not in product

                # Check for nitrile formation (forward synthesis)
                if not reactants_have_nitrile and product_has_nitrile:
                    has_nitrile_formation = True
                    print(f"Detected nitrile formation in reaction: {rsmi}")

                # Check for specific nitrile formation reactions
                if checker.check_reaction("Schmidt reaction nitrile", rsmi):
                    has_nitrile_formation = True
                    print(f"Detected Schmidt nitrile formation in reaction: {rsmi}")

                # Check for nitrile transformation (forward synthesis)
                if reactants_have_nitrile and not product_has_nitrile:
                    # Check for specific nitrile transformation reactions
                    if (
                        checker.check_reaction("Oxidation of nitrile to carboxylic acid", rsmi)
                        or checker.check_reaction("Nitrile to amide", rsmi)
                        or checker.check_reaction("Reduction of nitrile to amine", rsmi)
                        or checker.check_reaction("Grignard from nitrile to ketone", rsmi)
                    ):
                        has_nitrile_transformation = True
                        print(f"Detected specific nitrile transformation in reaction: {rsmi}")
                    else:
                        # Generic check for nitrile disappearance
                        has_nitrile_transformation = True
                        print(f"Detected potential nitrile transformation in reaction: {rsmi}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # If we found molecules with nitrile but didn't detect formation,
    # it suggests nitrile was formed at some point
    if molecules_with_nitrile and not has_nitrile_formation:
        # Check if nitrile appears in a non-starting material
        non_starting_nitriles = [
            m
            for m, d in molecules_with_nitrile
            if not any(
                node.get("in_stock", False)
                for node in route.get("children", [])
                if node["type"] == "mol" and node["smiles"] == m
            )
        ]
        if non_starting_nitriles:
            has_nitrile_formation = True
            print(
                f"Inferred nitrile formation from presence in non-starting materials: {non_starting_nitriles}"
            )

    # Return True if both formation and transformation of nitrile are detected
    strategy_present = has_nitrile_formation and has_nitrile_transformation
    print(f"Nitrile intermediate strategy detected: {strategy_present}")
    print(f"Nitrile formation detected: {has_nitrile_formation}")
    print(f"Nitrile transformation detected: {has_nitrile_transformation}")
    return strategy_present
