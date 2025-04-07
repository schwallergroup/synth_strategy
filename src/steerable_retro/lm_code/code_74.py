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
    This function detects if the synthetic route involves heterocycle construction,
    specifically focusing on benzofuran and other heterocycle formations.
    """
    heterocycle_construction_detected = False

    # List of heterocyclic rings to check for construction
    heterocycles = [
        "benzofuran",
        "benzothiophene",
        "indole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "furan",
        "pyrrole",
        "thiophene",
        "oxazole",
        "thiazole",
        "imidazole",
        "pyridine",
    ]

    # List of heterocycle formation reaction types
    heterocycle_reactions = [
        "{benzofuran}",
        "{benzothiophene}",
        "{indole}",
        "{benzoxazole}",
        "{benzothiazole}",
        "{benzimidazole_derivatives_aldehyde}",
        "{benzimidazole_derivatives_carboxylic-acid/ester}",
        "{thiazole}",
        "Paal-Knorr pyrrole synthesis",
        "{Paal-Knorr pyrrole}",
        "{Fischer indole}",
        "Formation of NOS Heterocycles",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_construction_detected

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if the reaction is a known heterocycle formation reaction
                for reaction_type in heterocycle_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Heterocycle formation reaction detected: {reaction_type}")
                        heterocycle_construction_detected = True
                        return

                # Check for heterocycle formation by comparing reactants and products
                for heterocycle in heterocycles:
                    # Check if product contains the heterocycle
                    product_has_heterocycle = checker.check_ring(heterocycle, product_smiles)

                    if product_has_heterocycle:
                        print(f"Product contains {heterocycle}: {product_smiles}")

                        # Check if any reactant contains the heterocycle
                        reactants_with_heterocycle = sum(
                            1
                            for r_smiles in reactants_smiles
                            if checker.check_ring(heterocycle, r_smiles)
                        )

                        # If product has heterocycle but reactants don't, it's a heterocycle construction
                        if reactants_with_heterocycle == 0:
                            print(
                                f"Heterocycle construction detected: {heterocycle} formed in reaction: {rsmi}"
                            )
                            heterocycle_construction_detected = True
                            return

                # Check for specific benzofuran formation patterns
                if checker.check_ring("benzofuran", product_smiles):
                    # Check for specific functional groups that might indicate benzofuran formation
                    if any(checker.check_fg("Phenol", r) for r in reactants_smiles) and any(
                        checker.check_fg("Primary halide", r)
                        or checker.check_fg("Secondary halide", r)
                        or checker.check_fg("Tertiary halide", r)
                        for r in reactants_smiles
                    ):
                        print(f"Benzofuran formation detected from phenol and halide: {rsmi}")
                        heterocycle_construction_detected = True
                        return
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting heterocycle construction detection...")
    dfs_traverse(route)
    print(f"Heterocycle construction detected: {heterocycle_construction_detected}")
    return heterocycle_construction_detected
