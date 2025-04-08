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
    This function detects if the core heterocyclic scaffold is formed early in the synthesis.
    """
    heterocycle_formation_detected = False
    max_depth = 0

    # List of common heterocycles to check
    heterocycles = [
        "pyridine",
        "pyrrole",
        "furan",
        "thiophene",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
        "morpholine",
        "piperidine",
        "piperazine",
        "oxirane",
        "aziridine",
        "thiirane",
    ]

    # First pass to determine the maximum depth
    def get_max_depth(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)
        for child in node.get("children", []):
            get_max_depth(child, depth + 1)

    get_max_depth(route)
    print(f"Maximum depth of synthesis tree: {max_depth}")

    # Second pass to detect heterocycle formation
    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_detected

        if node["type"] == "reaction" and depth >= max_depth / 2:  # Early in synthesis (high depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for heterocycle formation
                if product and Chem.MolFromSmiles(product):
                    product_mol = Chem.MolFromSmiles(product)
                    product_ring_count = len(product_mol.GetRingInfo().AtomRings())

                    # Check if product contains a heterocycle
                    product_has_heterocycle = False
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, product):
                            print(f"Found heterocycle {heterocycle} in product: {product}")
                            product_has_heterocycle = True
                            break

                    if product_has_heterocycle:
                        # Check if any reactant has fewer rings or doesn't have the heterocycle
                        for reactant in reactants:
                            if reactant and Chem.MolFromSmiles(reactant):
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                reactant_ring_count = len(reactant_mol.GetRingInfo().AtomRings())

                                # Check if reactant has the same heterocycle
                                reactant_has_heterocycle = False
                                for heterocycle in heterocycles:
                                    if checker.check_ring(heterocycle, reactant):
                                        reactant_has_heterocycle = True
                                        break

                                # If reactant has fewer rings or doesn't have the heterocycle,
                                # this indicates heterocycle formation
                                if (
                                    reactant_ring_count < product_ring_count
                                    or not reactant_has_heterocycle
                                ):
                                    print(f"Heterocycle formation detected at depth {depth}")
                                    print(f"Reactant: {reactant} (rings: {reactant_ring_count})")
                                    print(f"Product: {product} (rings: {product_ring_count})")
                                    heterocycle_formation_detected = True
                                    break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    return heterocycle_formation_detected
