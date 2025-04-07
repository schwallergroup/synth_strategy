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
    Detects linear synthesis strategy with heterocyclic scaffolds preserved
    throughout the synthesis.
    """
    # Track reaction steps and branching
    reaction_count = 0
    max_reactants_per_step = 0
    heterocycle_preserved = True

    # Define common heterocyclic rings to check
    heterocyclic_rings = [
        "indole",
        "quinoline",
        "isoquinoline",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "pyrrole",
        "furan",
        "thiophene",
        "imidazole",
        "oxazole",
        "thiazole",
        "triazole",
        "tetrazole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "piperidine",
        "morpholine",
        "piperazine",
        "pyrrolidine",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, max_reactants_per_step, heterocycle_preserved

        if node["type"] == "reaction":
            reaction_count += 1

            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Count reactants
            num_reactants = len(reactants_smiles)
            max_reactants_per_step = max(max_reactants_per_step, num_reactants)

            # Check if heterocycles are preserved across the reaction
            product_heterocycles = []
            for ring in heterocyclic_rings:
                if checker.check_ring(ring, product_smiles):
                    product_heterocycles.append(ring)

            reactant_heterocycles = []
            for reactant in reactants_smiles:
                for ring in heterocyclic_rings:
                    if checker.check_ring(ring, reactant):
                        reactant_heterocycles.append(ring)

            # If product has heterocycles but none of the reactants do, this is heterocycle formation
            # If reactants have heterocycles but product doesn't, this is heterocycle destruction
            if (product_heterocycles and not reactant_heterocycles) or (
                reactant_heterocycles and not product_heterocycles
            ):
                print(f"Heterocycle not preserved in reaction: {rsmi}")
                print(f"Product heterocycles: {product_heterocycles}")
                print(f"Reactant heterocycles: {reactant_heterocycles}")
                heterocycle_preserved = False

        elif node["type"] == "mol" and depth == 0:
            # Check if the final product contains any heterocycle
            final_product_smiles = node["smiles"]
            has_heterocycle = False
            for ring in heterocyclic_rings:
                if checker.check_ring(ring, final_product_smiles):
                    has_heterocycle = True
                    print(f"Final product contains heterocycle: {ring}")
                    break

            if not has_heterocycle:
                print("Final product does not contain any heterocycle")
                heterocycle_preserved = False

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Linear synthesis typically has max 2 reactants per step
    is_linear = max_reactants_per_step <= 2 and heterocycle_preserved and reaction_count > 0

    if is_linear:
        print(f"Found linear heterocycle synthesis with {reaction_count} steps")
    else:
        print(
            f"Not a linear heterocycle synthesis: max_reactants={max_reactants_per_step}, heterocycle_preserved={heterocycle_preserved}, reaction_count={reaction_count}"
        )

    return is_linear
