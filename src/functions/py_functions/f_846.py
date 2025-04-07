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
    This function detects a late-stage esterification in the synthetic route,
    specifically looking for conversion of carboxylic acid to ester in the final steps.
    """
    print("Starting analysis for late-stage esterification")
    esterification_depths = []

    def dfs_traverse(node, depth=0):
        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Analyzing reaction SMILES: {rsmi}")

                # Check if this is an esterification reaction
                is_esterification = checker.check_reaction(
                    "Esterification of Carboxylic Acids", rsmi
                ) or checker.check_reaction("Transesterification", rsmi)

                if is_esterification:
                    print(f"Found esterification reaction at depth {depth}")
                    esterification_depths.append(depth)
                else:
                    # Alternative check in case the reaction type is not directly recognized
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]
                    reactants_combined = ".".join(reactants)

                    # Check for both forward and retrosynthetic directions
                    # Forward: acid in reactants, ester in product
                    acid_in_reactants = any(
                        checker.check_fg("Carboxylic acid", reactant)
                        for reactant in reactants
                    )
                    ester_in_product = checker.check_fg("Ester", product)

                    # Retrosynthetic: ester in reactants, acid in product
                    ester_in_reactants = any(
                        checker.check_fg("Ester", reactant) for reactant in reactants
                    )
                    acid_in_product = checker.check_fg("Carboxylic acid", product)

                    if (acid_in_reactants and ester_in_product) or (
                        ester_in_reactants and acid_in_product
                    ):
                        print(
                            f"Found potential esterification at depth {depth} (acid ↔ ester conversion)"
                        )
                        esterification_depths.append(depth)
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if there's an esterification at depth 0, 1, or 2 (late stage)
    # Depth 0 is target molecule, depth 1 is first reaction, depth 2 is second reaction
    result = any(depth <= 2 for depth in esterification_depths)
    print(f"Esterification depths found: {esterification_depths}")
    print(f"Has late-stage esterification: {result}")
    return result
