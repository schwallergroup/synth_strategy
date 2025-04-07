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
    This function detects a linear synthesis strategy focused on fluorinated
    heterocycles without convergent steps.
    """
    # Track synthesis characteristics
    step_count = 0
    fluorinated_heterocycle_count = 0
    max_reactants_per_step = 0

    # List of common heterocycles to check
    heterocycles = [
        "pyridine",
        "pyrrole",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "triazole",
        "tetrazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "furan",
        "thiophene",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "piperidine",
        "piperazine",
        "morpholine",
        "pyrrolidine",
    ]

    def dfs_traverse(node):
        nonlocal step_count, fluorinated_heterocycle_count, max_reactants_per_step

        if node["type"] == "reaction":
            step_count += 1
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Count reactants in this step
            reactant_count = len(reactants)
            max_reactants_per_step = max(max_reactants_per_step, reactant_count)

            # Check for fluorinated heterocycles in the product
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                # Check if product contains a heterocycle
                contains_heterocycle = False
                for ring in heterocycles:
                    if checker.check_ring(ring, product):
                        contains_heterocycle = True
                        print(f"Found heterocycle: {ring} in product: {product}")
                        break

                # Check if product contains fluorine
                contains_fluorine = (
                    checker.check_fg("Trifluoro group", product)
                    or checker.check_fg("Aromatic halide", product)
                    and "[F]" in product
                    or checker.check_fg("Primary halide", product)
                    and "[F]" in product
                    or checker.check_fg("Secondary halide", product)
                    and "[F]" in product
                    or checker.check_fg("Tertiary halide", product)
                    and "[F]" in product
                    or checker.check_fg("Alkenyl halide", product)
                    and "[F]" in product
                    or checker.check_fg("Triflate", product)
                )

                if contains_heterocycle and contains_fluorine:
                    fluorinated_heterocycle_count += 1
                    print(f"Found fluorinated heterocycle (count: {fluorinated_heterocycle_count})")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Linear synthesis typically has few reactants per step and no convergent steps
    is_linear = max_reactants_per_step <= 2

    print(f"Step count: {step_count}, Max reactants per step: {max_reactants_per_step}")
    print(f"Fluorinated heterocycle count: {fluorinated_heterocycle_count}")
    print(f"Is linear synthesis: {is_linear}")

    # Return true if it's a linear synthesis with at least one fluorinated heterocycle
    return is_linear and fluorinated_heterocycle_count > 0
