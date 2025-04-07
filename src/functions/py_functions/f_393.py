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
    Detects a strategy involving the formation of a heterocyclic system,
    specifically looking for pyrazolopyrimidine scaffold construction.
    """
    # Track if we found heterocycle formation
    heterocycle_formed = False

    # List of heterocyclic rings to check
    heterocyclic_rings = [
        "pyrazole",
        "pyrimidine",
        "triazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "isoxazole",
        "isothiazole",
        "tetrazole",
        "furan",
        "thiophene",
        "pyrrole",
        "pyridine",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formed

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}")
                print(f"Product: {product_smiles}")
                print(f"Reactants: {reactants_smiles}")

                # Check for heterocycle formation
                for ring in heterocyclic_rings:
                    # Check if product contains the heterocycle
                    if checker.check_ring(ring, product_smiles):
                        print(f"Found {ring} in product")

                        # Check if any reactant contains the same heterocycle
                        reactant_has_ring = any(
                            checker.check_ring(ring, r) for r in reactants_smiles
                        )

                        if not reactant_has_ring:
                            print(
                                f"Heterocycle formation detected: {ring} at depth {depth}"
                            )
                            heterocycle_formed = True

                            # Check specifically for pyrazolopyrimidine-like structures
                            if ring in ["pyrazole", "pyrimidine"] and any(
                                checker.check_reaction(rxn, rsmi)
                                for rxn in [
                                    "Formation of NOS Heterocycles",
                                    "{pyrazole}",
                                    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                                    "Huisgen 1,3 dipolar cycloaddition",
                                ]
                            ):
                                print(
                                    f"Confirmed heterocycle formation reaction at depth {depth}"
                                )
                                return  # Found what we're looking for

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if heterocycle_formed:
        print("Heterocycle formation strategy detected")
    else:
        print("No heterocycle formation strategy detected")

    return heterocycle_formed
