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
    Detects the functional group transformation sequence: ketone -> oxime -> isoxazole
    """
    # Dictionary to track molecules through transformations
    # Key: molecule SMILES, Value: transformation stage (0=ketone, 1=oxime, 2=isoxazole)
    transformation_stages = {}
    sequence_completed = False

    def dfs_traverse(node, depth=0):
        nonlocal transformation_stages, sequence_completed

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # In retrosynthesis, product is the starting material and reactants are the targets
                # Check if product contains isoxazole
                if checker.check_ring("isoxazole", product):
                    print(f"Found isoxazole in product: {product}")

                    # Check if any reactant contains oxime
                    for reactant in reactants:
                        if checker.check_fg("Oxime", reactant):
                            print(f"Found oxime in reactant: {reactant}")

                            # This is a transformation from oxime to isoxazole
                            transformation_stages[reactant] = 1  # oxime stage
                            transformation_stages[product] = 2  # isoxazole stage

                # Check if product contains oxime
                if checker.check_fg("Oxime", product):
                    print(f"Found oxime in product: {product}")

                    # Check if any reactant contains ketone
                    for reactant in reactants:
                        if checker.check_fg("Ketone", reactant):
                            print(f"Found ketone in reactant: {reactant}")

                            # This is a transformation from ketone to oxime
                            transformation_stages[reactant] = 0  # ketone stage
                            transformation_stages[product] = 1  # oxime stage

                            # Check if this oxime was later transformed to isoxazole
                            if (
                                product in transformation_stages
                                and transformation_stages[product] == 1
                            ):
                                for other_product, stage in transformation_stages.items():
                                    if stage == 2 and other_product != product:
                                        # We found the complete sequence
                                        print(
                                            f"Complete sequence found: {reactant} (ketone) -> {product} (oxime) -> {other_product} (isoxazole)"
                                        )
                                        sequence_completed = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return sequence_completed
