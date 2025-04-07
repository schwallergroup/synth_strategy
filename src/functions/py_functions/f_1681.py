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
    This function detects if the synthetic route follows a convergent approach
    where multiple heterocyclic fragments are combined rather than a linear synthesis.
    """
    fragment_combinations = 0

    # List of common heterocycles to check
    heterocycle_types = [
        "pyridine",
        "pyrrole",
        "furan",
        "thiophene",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrazole",
        "isoxazole",
        "isothiazole",
        "triazole",
        "tetrazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "indole",
        "benzofuran",
        "benzothiophene",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
        "piperidine",
        "tetrahydrofuran",
        "tetrahydropyran",
        "morpholine",
        "thiomorpholine",
    ]

    def dfs_traverse(node):
        nonlocal fragment_combinations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if multiple fragments are being combined
                if len(reactants) >= 2:
                    # Count heterocycles in each reactant
                    heterocycle_counts = []
                    for reactant in reactants:
                        heterocycle_count = 0
                        for heterocycle in heterocycle_types:
                            if checker.check_ring(heterocycle, reactant):
                                heterocycle_count += 1
                                break  # Count each reactant only once
                        heterocycle_counts.append(heterocycle_count)

                    # Check if product contains heterocycles
                    product_heterocycles = sum(
                        1
                        for heterocycle in heterocycle_types
                        if checker.check_ring(heterocycle, product)
                    )

                    # If at least 2 reactants have heterocycles and product has heterocycles
                    if (
                        sum(1 for count in heterocycle_counts if count > 0) >= 2
                        and product_heterocycles > 0
                    ):
                        fragment_combinations += 1
                        print(f"Heterocycle fragment combination detected: {rsmi}")
                        print(f"Heterocycle counts in reactants: {heterocycle_counts}")
                        print(f"Heterocycles in product: {product_heterocycles}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    is_convergent = fragment_combinations >= 1
    print(f"Fragment combinations: {fragment_combinations}")
    print(f"Convergent heterocycle synthesis: {is_convergent}")

    return is_convergent
