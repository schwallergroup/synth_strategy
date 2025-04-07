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
    This function detects convergent synthesis with heteroaromatic building blocks.
    """
    # Track if we found the pattern
    found_convergent_step = False
    heteroaromatic_fragments = 0

    # List of common heteroaromatic rings to check
    heteroaromatic_rings = [
        "pyridine",
        "pyrrole",
        "furan",
        "thiophene",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
        "purine",
        "triazole",
        "tetrazole",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_convergent_step, heteroaromatic_fragments

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for convergent step (multiple complex fragments)
                if len(reactants_smiles) >= 2:
                    # Count heteroaromatic fragments
                    heteroaromatic_count = 0
                    for reactant in reactants_smiles:
                        # Check if molecule has heteroaromatic rings
                        has_heteroaromatic_ring = False
                        for ring_name in heteroaromatic_rings:
                            if checker.check_ring(ring_name, reactant):
                                print(
                                    f"Found heteroaromatic ring {ring_name} in reactant {reactant}"
                                )
                                has_heteroaromatic_ring = True
                                break

                        if has_heteroaromatic_ring:
                            heteroaromatic_count += 1

                    if heteroaromatic_count >= 2:
                        found_convergent_step = True
                        heteroaromatic_fragments = heteroaromatic_count
                        print(
                            f"Found convergent step with {heteroaromatic_count} heteroaromatic fragments at depth {depth}"
                        )
                        print(f"Reaction SMILES: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found a convergent step with at least 2 heteroaromatic fragments
    return found_convergent_step and heteroaromatic_fragments >= 2
