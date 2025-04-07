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
    This function detects a strategy involving heterocycle formation:
    1. Formation of a heterocyclic ring system
    2. Occurs in the middle of the synthesis (not at the beginning or end)
    """

    # Track if we found heterocycle formation
    found_heterocycle_formation = False
    max_depth = 0

    # First pass to determine the maximum depth of the synthesis route
    def get_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            get_max_depth(child, current_depth + 1)

    get_max_depth(route)
    print(f"Maximum depth of synthesis route: {max_depth}")

    # Define heterocyclic rings to check for
    heterocyclic_rings = [
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazine",
        "furan",
        "thiophene",
        "pyrrole",
        "imidazole",
        "oxazole",
        "thiazole",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
        "quinazoline",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "pyrrolidine",
        "oxazolidine",
        "thiazolidine",
    ]

    def is_heterocycle_formation(reaction_node):
        """Check if a reaction forms a heterocycle"""
        try:
            rsmi = reaction_node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return False

            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants = reactants_part.split(".")
            product = product_part

            # Check if product contains any heterocyclic ring that's not in reactants
            product_mol = Chem.MolFromSmiles(product)
            if not product_mol:
                return False

            # Check for heterocycle formation
            for ring in heterocyclic_rings:
                # Check if product has the heterocycle
                if checker.check_ring(ring, product):
                    # Check if any reactant already has the heterocycle
                    reactant_has_ring = False
                    for reactant in reactants:
                        if checker.check_ring(ring, reactant):
                            reactant_has_ring = True
                            break

                    # If product has the ring but reactants don't, it's a formation
                    if not reactant_has_ring:
                        print(f"Found formation of {ring} heterocycle")
                        return True

            return False
        except Exception as e:
            print(f"Error in is_heterocycle_formation: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal found_heterocycle_formation

        if node["type"] == "reaction":
            # Define middle stage as between 25% and 75% of max depth
            lower_bound = max(1, int(max_depth * 0.25))
            upper_bound = int(max_depth * 0.75)

            # Check for heterocycle formation in middle of synthesis
            if lower_bound <= depth <= upper_bound and is_heterocycle_formation(node):
                print(f"Found heterocycle formation at depth {depth} (middle stage)")
                found_heterocycle_formation = True

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Heterocycle formation strategy detected: {found_heterocycle_formation}")
    return found_heterocycle_formation
