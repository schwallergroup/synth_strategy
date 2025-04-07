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
    This function detects a strategy involving the use of azide as a synthetic intermediate.
    """
    azide_formation = False
    azide_conversion = False

    def dfs_traverse(node):
        nonlocal azide_formation, azide_conversion

        if node["type"] == "reaction":
            metadata = node.get("metadata", {})
            rsmi = metadata.get("rsmi", "")
            if not rsmi:
                return

            # Split reaction SMILES into reactants and products
            parts = rsmi.split(">")
            if len(parts) >= 3:
                reactants_smiles = parts[0]
                product_smiles = parts[-1]

                # Check for azide formation
                if checker.check_fg("Azide", product_smiles) and not checker.check_fg(
                    "Azide", reactants_smiles
                ):
                    print(f"Found azide formation in reaction: {rsmi}")
                    azide_formation = True

                    # Check for specific azide formation reactions
                    if (
                        checker.check_reaction(
                            "Formation of Azides from halogens", rsmi
                        )
                        or checker.check_reaction(
                            "Formation of Azides from boronic acids", rsmi
                        )
                        or checker.check_reaction("Alcohol to azide", rsmi)
                        or checker.check_reaction("Amine to azide", rsmi)
                    ):
                        print(f"Confirmed azide formation reaction: {rsmi}")

                # Check for azide conversion
                if checker.check_fg("Azide", reactants_smiles) and not checker.check_fg(
                    "Azide", product_smiles
                ):
                    print(f"Found azide conversion in reaction: {rsmi}")
                    azide_conversion = True

                    # Check for specific azide conversion reactions
                    if (
                        checker.check_reaction(
                            "Azide to amine reduction (Staudinger)", rsmi
                        )
                        or checker.check_reaction(
                            "Huisgen alkyne-azide 1,3 dipolar cycloaddition", rsmi
                        )
                        or checker.check_reaction(
                            "Huisgen 1,3 dipolar cycloaddition", rsmi
                        )
                        or checker.check_reaction(
                            "Huisgen alkene-azide 1,3 dipolar cycloaddition", rsmi
                        )
                        or checker.check_reaction(
                            "Huisgen_Cu-catalyzed_1,4-subst", rsmi
                        )
                        or checker.check_reaction(
                            "Huisgen_Ru-catalyzed_1,5_subst", rsmi
                        )
                        or checker.check_reaction("Huisgen_disubst-alkyne", rsmi)
                    ):
                        print(f"Confirmed azide conversion reaction: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Azide formation detected: {azide_formation}")
    print(f"Azide conversion detected: {azide_conversion}")

    return azide_formation and azide_conversion
