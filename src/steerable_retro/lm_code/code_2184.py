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
    This function detects tetrazole ring formation via azide-isocyanate cycloaddition.
    """
    tetrazole_found = False
    azide_isocyanate_reaction = False

    def dfs_traverse(node):
        nonlocal tetrazole_found, azide_isocyanate_reaction

        if node["type"] == "reaction":
            # Check if this reaction forms a tetrazole
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains tetrazole
                if checker.check_ring("tetrazole", product):
                    print(f"Found tetrazole in product: {product}")
                    tetrazole_found = True

                    # Check if reactants contain azide and isocyanate
                    has_azide = False
                    has_isocyanate = False

                    for reactant in reactants:
                        if checker.check_fg("Azide", reactant):
                            print(f"Found azide in reactant: {reactant}")
                            has_azide = True
                        if checker.check_fg("Isocyanate", reactant):
                            print(f"Found isocyanate in reactant: {reactant}")
                            has_isocyanate = True

                    # Check if this is a Huisgen cycloaddition reaction
                    if (
                        checker.check_reaction("Huisgen 1,3 dipolar cycloaddition", rsmi)
                        or checker.check_reaction(
                            "Huisgen alkyne-azide 1,3 dipolar cycloaddition", rsmi
                        )
                        or checker.check_reaction(
                            "Azide-nitrile click cycloaddition to tetrazole", rsmi
                        )
                    ):
                        print(f"Confirmed cycloaddition reaction: {rsmi}")
                        azide_isocyanate_reaction = True
                    elif has_azide and has_isocyanate:
                        print(
                            f"Found azide and isocyanate in reactants for tetrazole formation: {rsmi}"
                        )
                        azide_isocyanate_reaction = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(
        f"tetrazole_found: {tetrazole_found}, azide_isocyanate_reaction: {azide_isocyanate_reaction}"
    )
    return tetrazole_found and azide_isocyanate_reaction
