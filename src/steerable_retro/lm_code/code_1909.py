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
    Detects if the synthetic route involves late-stage tetrazole formation from nitrile
    via [3+2] cycloaddition with azide.
    """
    tetrazole_formed = False
    nitrile_precursor = False
    azide_precursor = False
    tetrazole_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal tetrazole_formed, nitrile_precursor, azide_precursor, tetrazole_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this reaction forms a tetrazole via azide-nitrile cycloaddition
            if checker.check_reaction("Azide-nitrile click cycloaddition to tetrazole", rsmi):
                print(f"Found tetrazole formation reaction at depth {depth}: {rsmi}")
                tetrazole_formed = True
                tetrazole_depth = min(tetrazole_depth, depth)

                # Check if reactants contain nitrile and azide
                for reactant in reactants:
                    if checker.check_fg("Nitrile", reactant):
                        print(f"Found nitrile precursor: {reactant}")
                        nitrile_precursor = True
                    if checker.check_fg("Azide", reactant):
                        print(f"Found azide precursor: {reactant}")
                        azide_precursor = True

            # Alternative check: product has tetrazole and reactants have nitrile/azide
            elif checker.check_ring("tetrazole", product):
                print(f"Found tetrazole in product at depth {depth}: {product}")
                has_nitrile = False
                has_azide = False

                for reactant in reactants:
                    if checker.check_fg("Nitrile", reactant):
                        print(f"Found nitrile precursor: {reactant}")
                        has_nitrile = True
                        nitrile_precursor = True
                    if checker.check_fg("Azide", reactant):
                        print(f"Found azide precursor: {reactant}")
                        has_azide = True
                        azide_precursor = True

                # If we have both nitrile and azide, it's definitely tetrazole formation
                if has_nitrile and has_azide:
                    tetrazole_formed = True
                    tetrazole_depth = min(tetrazole_depth, depth)
                # If we only have nitrile but product has tetrazole, it's likely tetrazole formation
                elif has_nitrile:
                    tetrazole_formed = True
                    tetrazole_depth = min(tetrazole_depth, depth)
                    print(
                        "Inferred tetrazole formation from nitrile (azide not explicitly detected)"
                    )

        # For molecule nodes, check if they contain tetrazole
        elif node["type"] == "mol":
            if checker.check_ring("tetrazole", node["smiles"]):
                print(f"Found tetrazole in molecule at depth {depth}: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if tetrazole formation is late-stage (depth <= 2)
    is_late_stage = tetrazole_depth <= 2

    print(f"Tetrazole formed: {tetrazole_formed}")
    print(f"From nitrile: {nitrile_precursor}")
    print(f"With azide: {azide_precursor}")
    print(f"Tetrazole formation depth: {tetrazole_depth}")
    print(f"Is late-stage: {is_late_stage}")

    # For this specific function, we'll consider it valid if we have tetrazole formation
    # from nitrile at a late stage, even if azide isn't explicitly detected
    return tetrazole_formed and nitrile_precursor and is_late_stage
