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
    This function detects if the synthetic route involves late-stage heterocycle formation,
    specifically a pyrazole ring formation in the second half of the synthesis.
    """
    heterocycle_formation_found = False
    heterocycle_formation_depth = -1
    max_depth = -1

    # List of heterocycles to check
    heterocycles = [
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "isoxazole",
        "isothiazole",
        "triazole",
        "tetrazole",
        "oxadiazole",
        "thiadiazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_found, heterocycle_formation_depth, max_depth

        # Update max depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check for heterocycle formation
            product_has_heterocycle = False
            reactants_have_heterocycle = False
            formed_heterocycle = None

            # Check which heterocycle is formed
            for heterocycle in heterocycles:
                if checker.check_ring(heterocycle, product_smiles):
                    print(f"Product contains {heterocycle} ring")
                    product_has_heterocycle = True
                    formed_heterocycle = heterocycle

                    # Check if any reactant has the same heterocycle
                    for reactant in reactants_smiles:
                        if checker.check_ring(heterocycle, reactant):
                            print(f"Reactant also contains {heterocycle} ring")
                            reactants_have_heterocycle = True
                            break

                    # If heterocycle is in product but not in reactants, it's formed in this reaction
                    if not reactants_have_heterocycle:
                        # Check if this is a known heterocycle formation reaction
                        reaction_is_heterocycle_formation = False

                        # Check for specific reaction types that form heterocycles
                        reaction_names_to_check = [
                            "pyrazole formation",
                            "pyrazole",
                            "{pyrazole}",
                            "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                            "Huisgen 1,3 dipolar cycloaddition",
                            "Huisgen alkene-azide 1,3 dipolar cycloaddition",
                            "[3+2]-cycloaddition of hydrazone and alkyne",
                            "[3+2]-cycloaddition of hydrazone and alkene",
                            "[3+2]-cycloaddition of diazoalkane and alkyne",
                            "[3+2]-cycloaddition of diazoalkane and alkene",
                        ]

                        for reaction_name in reaction_names_to_check:
                            if checker.check_reaction(reaction_name, rsmi):
                                print(f"Confirmed {reaction_name} formation reaction")
                                reaction_is_heterocycle_formation = True
                                break

                        # Direct check for pyrazole formation
                        if heterocycle == "pyrazole" and not reaction_is_heterocycle_formation:
                            # Check for hydrazine in reactants and nitrile in reactants
                            has_hydrazine = False
                            has_nitrile = False

                            for reactant in reactants_smiles:
                                if checker.check_fg("Hydrazine", reactant):
                                    print("Found hydrazine in reactants")
                                    has_hydrazine = True
                                if checker.check_fg("Nitrile", reactant):
                                    print("Found nitrile in reactants")
                                    has_nitrile = True

                            if has_hydrazine and has_nitrile:
                                print("Detected pyrazole formation from hydrazine and nitrile")
                                reaction_is_heterocycle_formation = True

                        if reaction_is_heterocycle_formation:
                            print(f"{formed_heterocycle} ring formation detected at depth {depth}")
                            heterocycle_formation_found = True
                            # If multiple formations are found, keep the one with the lowest depth (latest stage)
                            if (
                                heterocycle_formation_depth == -1
                                or depth < heterocycle_formation_depth
                            ):
                                heterocycle_formation_depth = depth

                    break  # Once we find a heterocycle in the product, no need to check others

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Max depth: {max_depth}")
    print(f"Heterocycle formation found: {heterocycle_formation_found}")
    print(f"Heterocycle formation depth: {heterocycle_formation_depth}")

    # Calculate the halfway point of the synthesis
    halfway_depth = max_depth / 2
    print(f"Halfway depth: {halfway_depth}")

    # Check if heterocycle formation occurred in the second half of synthesis
    # In retrosynthetic analysis, lower depth values correspond to later stages in synthesis
    # So we need to check if the formation depth is less than halfway_depth
    is_late_stage = heterocycle_formation_found and heterocycle_formation_depth < halfway_depth
    print(f"Is late stage heterocycle formation: {is_late_stage}")

    return is_late_stage
