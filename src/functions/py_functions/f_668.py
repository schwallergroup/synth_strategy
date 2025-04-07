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
    Detects a synthetic strategy involving the formation of a sulfonamide
    from a sulfonyl chloride and an amine.
    """
    found_sulfonamide_formation = False

    def dfs_traverse(node):
        nonlocal found_sulfonamide_formation

        if node["type"] == "reaction" and not found_sulfonamide_formation:
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # First check if this is a known sulfonamide synthesis reaction
                if checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                ) or checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                ):
                    found_sulfonamide_formation = True
                    print(
                        "Found sulfonamide formation step (Schotten-Baumann reaction)"
                    )
                    return

                # If not a known reaction type, check for the functional group transformation
                sulfonyl_halide_in_reactants = any(
                    checker.check_fg("Sulfonyl halide", reactant)
                    for reactant in reactants
                )
                sulfonamide_in_product = checker.check_fg("Sulfonamide", product)

                if sulfonyl_halide_in_reactants and sulfonamide_in_product:
                    # Check if any reactant has an amine
                    for reactant in reactants:
                        if checker.check_fg(
                            "Primary amine", reactant
                        ) or checker.check_fg("Secondary amine", reactant):
                            found_sulfonamide_formation = True
                            print(
                                "Found sulfonamide formation step (functional group transformation)"
                            )
                            return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_sulfonamide_formation
