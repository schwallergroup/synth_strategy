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
    This function detects heterocyclic ring formation, particularly pyrazole formation
    via cyclization of a hydrazine derivative.
    """
    ring_formation_found = False

    def dfs_traverse(node):
        nonlocal ring_formation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Check if this is a pyrazole formation reaction
                if checker.check_reaction("pyrazole", rsmi):
                    print(f"Found pyrazole formation reaction: {rsmi}")
                    ring_formation_found = True
                    return

                # If not directly identified as pyrazole formation, check for hydrazine cyclization
                reactants = reactants_str.split(".")

                # Check if any reactant contains hydrazine group
                hydrazine_present = any(checker.check_fg("Hydrazine", r) for r in reactants)

                # Check if product contains pyrazole ring
                pyrazole_present = checker.check_ring("pyrazole", product_str)

                if hydrazine_present and pyrazole_present:
                    print(f"Found pyrazole ring formation via hydrazine cyclization: {rsmi}")
                    ring_formation_found = True
                    return

                # Check for other heterocyclic ring formations
                # Check if product contains a heterocyclic ring that wasn't in any reactant
                heterocyclic_rings = [
                    "pyrazole",
                    "imidazole",
                    "oxazole",
                    "thiazole",
                    "triazole",
                    "tetrazole",
                ]

                for ring in heterocyclic_rings:
                    if checker.check_ring(ring, product_str):
                        # Check if the ring was present in any reactant
                        ring_in_reactants = any(checker.check_ring(ring, r) for r in reactants)

                        if not ring_in_reactants:
                            print(f"Found {ring} ring formation: {rsmi}")
                            ring_formation_found = True
                            return

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return ring_formation_found
