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
    This function detects if the synthesis includes a mid-stage heterocycle formation.
    """
    has_heterocycle_formation = False
    total_depth = 0

    # List of common heterocycles to check
    heterocycles = [
        "imidazole",
        "pyrrole",
        "pyrazole",
        "oxazole",
        "thiazole",
        "triazole",
        "tetrazole",
        "furan",
        "thiophene",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "indazole",
        "benzotriazole",
    ]

    # First pass to determine total depth
    def get_max_depth(node, current_depth=0):
        nonlocal total_depth
        total_depth = max(total_depth, current_depth)

        for child in node.get("children", []):
            get_max_depth(child, current_depth + 1)

    # Second pass to check for heterocycle formation
    def dfs_traverse(node, depth=0):
        nonlocal has_heterocycle_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # We consider "mid-stage" to be in the middle third of the synthesis
            is_mid_stage = depth >= total_depth // 3 and depth <= (2 * total_depth) // 3

            if is_mid_stage:
                try:
                    # Extract reactants and product
                    rsmi = node["metadata"]["rsmi"]
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    print(f"Checking reaction at depth {depth}/{total_depth}: {rsmi}")

                    # Check for heterocycle formation
                    for heterocycle in heterocycles:
                        # Check if heterocycle is in product
                        if checker.check_ring(heterocycle, product_smiles):
                            print(f"Found {heterocycle} in product: {product_smiles}")

                            # Check if heterocycle wasn't present in any reactant
                            heterocycle_in_reactants = False
                            for reactant_smiles in reactants_smiles:
                                if checker.check_ring(heterocycle, reactant_smiles):
                                    heterocycle_in_reactants = True
                                    print(f"Found {heterocycle} in reactant: {reactant_smiles}")
                                    break

                            if not heterocycle_in_reactants:
                                print(
                                    f"Detected mid-stage {heterocycle} formation at depth {depth}"
                                )
                                has_heterocycle_formation = True
                                break
                except Exception as e:
                    print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Get total depth first
    get_max_depth(route)
    print(f"Total synthesis depth: {total_depth}")

    # Start traversal
    dfs_traverse(route)

    return has_heterocycle_formation
