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
    This function detects the use of azide chemistry as a key strategy in the synthesis,
    specifically looking for azide formation and utilization.
    """
    # Track azide-related transformations
    azide_forming_reactions = 0
    azide_utilizing_reactions = 0
    azide_intermediates = 0

    def dfs_traverse(node):
        nonlocal azide_forming_reactions, azide_utilizing_reactions, azide_intermediates

        if node["type"] == "mol":
            # Check if molecule contains azide group
            if checker.check_fg("Azide", node["smiles"]):
                azide_intermediates += 1
                print(f"Azide intermediate detected: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for azide formation reactions
            if (
                checker.check_reaction("Formation of Azides from halogens", rsmi)
                or checker.check_reaction("Formation of Azides from boronic acids", rsmi)
                or checker.check_reaction("Amine to azide", rsmi)
            ):
                azide_forming_reactions += 1
                print(f"Azide formation reaction detected: {rsmi}")

            # Check for azide utilization reactions
            elif (
                checker.check_reaction("Huisgen alkyne-azide 1,3 dipolar cycloaddition", rsmi)
                or checker.check_reaction("Huisgen 1,3 dipolar cycloaddition", rsmi)
                or checker.check_reaction("Huisgen alkene-azide 1,3 dipolar cycloaddition", rsmi)
                or checker.check_reaction("Azide to amine reduction (Staudinger)", rsmi)
                or checker.check_reaction("Azide-nitrile click cycloaddition to tetrazole", rsmi)
                or checker.check_reaction("Azide-nitrile click cycloaddition to triazole", rsmi)
            ):
                azide_utilizing_reactions += 1
                print(f"Azide utilizing reaction detected: {rsmi}")

            # Fallback method if specific reaction checks fail
            else:
                # Check for azide formation
                reactants_have_precursor = any(
                    checker.check_fg("Primary amine", r)
                    or checker.check_fg("Secondary amine", r)
                    or checker.check_fg("Tertiary amine", r)
                    or checker.check_fg("Primary halide", r)
                    or checker.check_fg("Secondary halide", r)
                    or checker.check_fg("Tertiary halide", r)
                    or checker.check_fg("Aromatic halide", r)
                    or checker.check_fg("Boronic acid", r)
                    or checker.check_fg("Boronic ester", r)
                    for r in reactants
                )

                product_has_azide = checker.check_fg("Azide", product)

                if reactants_have_precursor and product_has_azide:
                    azide_forming_reactions += 1
                    print(f"Azide formation detected through functional group analysis: {rsmi}")

                # Check for azide utilization
                reactants_have_azide = any(checker.check_fg("Azide", r) for r in reactants)
                product_lacks_azide = not checker.check_fg("Azide", product)

                if reactants_have_azide and product_lacks_azide:
                    azide_utilizing_reactions += 1
                    print(f"Azide utilization detected through functional group analysis: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Determine if our strategy is present
    # Strategy is present if we have azide intermediates AND either formation or utilization reactions
    azide_strategy_present = azide_intermediates > 0 and (
        azide_forming_reactions > 0 or azide_utilizing_reactions > 0
    )

    print(f"Azide-based functional group manipulation strategy detected: {azide_strategy_present}")
    print(f"Azide forming reactions: {azide_forming_reactions}")
    print(f"Azide utilizing reactions: {azide_utilizing_reactions}")
    print(f"Azide intermediates: {azide_intermediates}")

    return azide_strategy_present
