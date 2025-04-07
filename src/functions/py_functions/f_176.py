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
    Detects if the synthetic route involves building complexity on a
    pyridine-containing scaffold.
    """
    # Track if we found a reaction that builds a heterocycle on a pyridine scaffold
    heterocycle_elaboration_found = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_elaboration_found

        # Check reaction nodes
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if reactants contain pyridine
            reactant_has_pyridine = False
            for reactant in reactants_part.split("."):
                if checker.check_ring("pyridine", reactant):
                    reactant_has_pyridine = True
                    break

            # Check if product has a more complex heterocycle structure
            product_has_complex_heterocycle = False
            product_has_pyridine = checker.check_ring("pyridine", product_part)

            # Check for various heterocycles that could be built on pyridine
            if product_has_pyridine:
                # Check for imidazopyridine-like structures (pyridine fused with other heterocycles)
                if (
                    checker.check_ring("imidazole", product_part)
                    or checker.check_ring("oxazole", product_part)
                    or checker.check_ring("thiazole", product_part)
                    or checker.check_ring("pyrazole", product_part)
                    or checker.check_ring("triazole", product_part)
                ):
                    product_has_complex_heterocycle = True

                # Check for quinoline or isoquinoline (pyridine fused with benzene)
                if checker.check_ring("quinoline", product_part) or checker.check_ring(
                    "isoquinoline", product_part
                ):
                    product_has_complex_heterocycle = True

            # If the reaction builds a complex heterocycle on a pyridine scaffold
            if reactant_has_pyridine and product_has_complex_heterocycle:
                print(f"Found heterocycle elaboration reaction: {rsmi}")
                heterocycle_elaboration_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return heterocycle_elaboration_found
