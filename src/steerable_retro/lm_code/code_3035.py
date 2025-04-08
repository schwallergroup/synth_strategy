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
    Detects if the synthesis route involves late-stage imidazole ring formation.
    Late stage means at low depth in the tree (closer to the final product).
    """
    imidazole_formed = False
    imidazole_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal imidazole_formed, imidazole_depth

        if node["type"] == "reaction":
            # Extract product and reactants
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if imidazole is formed in this reaction
            product_has_imidazole = checker.check_ring("imidazole", product)

            reactants_have_imidazole = False
            for reactant in reactants:
                if checker.check_ring("imidazole", reactant):
                    reactants_have_imidazole = True
                    break

            # Check if this reaction forms an imidazole ring
            if product_has_imidazole and not reactants_have_imidazole:
                # Check if this is a known imidazole-forming reaction
                imidazole_reaction = (
                    checker.check_reaction("triaryl-imidazole", rsmi)
                    or checker.check_reaction("imidazole", rsmi)
                    or checker.check_reaction("{imidazole}", rsmi)
                )

                if imidazole_reaction:
                    imidazole_formed = True
                    imidazole_depth = min(imidazole_depth, depth)
                    print(f"Imidazole formation detected at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")
                    print(f"Product: {product}")
                    print(f"Reactants: {reactants}")
                else:
                    # Even if not a known reaction type, still consider it if imidazole appears
                    imidazole_formed = True
                    imidazole_depth = min(imidazole_depth, depth)
                    print(f"Imidazole formation detected (unknown reaction type) at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Consider it late-stage if it happens at depth 0, 1, or 2
    return imidazole_formed and imidazole_depth <= 2
