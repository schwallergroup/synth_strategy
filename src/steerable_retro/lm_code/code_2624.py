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
    This function detects a strategy involving transformation of a nitro group
    during heterocycle formation.
    """
    nitro_in_heterocycle_formation = False

    # List of common heterocyclic rings to check
    heterocycles = [
        "furan",
        "pyrrole",
        "thiophene",
        "pyridine",
        "oxazole",
        "thiazole",
        "imidazole",
        "pyrazole",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "triazole",
        "tetrazole",
        "pyrimidine",
        "pyrazine",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "indole",
        "quinoline",
        "isoquinoline",
    ]

    def dfs_traverse(node):
        nonlocal nitro_in_heterocycle_formation

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for nitro group in reactants
                has_nitro_in_reactants = any(
                    checker.check_fg("Nitro group", r) for r in reactants_smiles
                )

                if has_nitro_in_reactants:
                    print(f"Found nitro group in reactants: {rsmi}")

                    # Check if nitro group is absent in product (transformed)
                    nitro_transformed = not checker.check_fg("Nitro group", product_smiles)

                    if nitro_transformed:
                        print(f"Nitro group transformed in product: {product_smiles}")

                        # Check for heterocycle formation reaction
                        is_heterocycle_formation = checker.check_reaction(
                            "Formation of NOS Heterocycles", rsmi
                        )

                        # If not a direct heterocycle formation reaction, check for new heterocycles in product
                        if not is_heterocycle_formation:
                            # Check for heterocycles in product
                            product_heterocycles = [
                                cycle
                                for cycle in heterocycles
                                if checker.check_ring(cycle, product_smiles)
                            ]

                            # Check for heterocycles in reactants
                            reactant_heterocycles = []
                            for r in reactants_smiles:
                                reactant_heterocycles.extend(
                                    [
                                        cycle
                                        for cycle in heterocycles
                                        if checker.check_ring(cycle, r)
                                    ]
                                )

                            # Check if new heterocycles formed
                            new_heterocycles = [
                                cycle
                                for cycle in product_heterocycles
                                if cycle not in reactant_heterocycles
                            ]

                            if new_heterocycles:
                                print(f"New heterocycles formed: {new_heterocycles}")
                                is_heterocycle_formation = True

                        if is_heterocycle_formation:
                            print(
                                f"Detected nitro group transformation in heterocycle formation: {rsmi}"
                            )
                            nitro_in_heterocycle_formation = True
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return nitro_in_heterocycle_formation
