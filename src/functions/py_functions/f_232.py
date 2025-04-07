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
    This function detects the formation of a quinazoline core during the synthesis.
    """
    from rdkit import Chem

    quinazoline_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal quinazoline_formation_detected

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a Niementowski quinazoline reaction
                if checker.check_reaction("{Niementowski_quinazoline}", rsmi):
                    print(
                        f"Detected Niementowski quinazoline reaction at depth {depth}: {rsmi}"
                    )
                    quinazoline_formation_detected = True
                    return

                # Check if product contains a quinazoline structure
                try:
                    product_has_quinazoline = checker.check_ring(
                        "quinazoline", product_smiles
                    )

                    if product_has_quinazoline:
                        print(
                            f"Product contains quinazoline at depth {depth}: {product_smiles}"
                        )

                        # Check if any reactant has quinazoline
                        reactant_has_quinazoline = False
                        for reactant_smiles in reactants_smiles:
                            if checker.check_ring("quinazoline", reactant_smiles):
                                reactant_has_quinazoline = True
                                print(
                                    f"Reactant already contains quinazoline: {reactant_smiles}"
                                )
                                break

                        if not reactant_has_quinazoline:
                            print(
                                f"Quinazoline core formation detected in reaction at depth {depth}: {rsmi}"
                            )
                            quinazoline_formation_detected = True
                except Exception as e:
                    print(f"Error checking quinazoline structure: {e}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)
    print(f"Quinazoline formation detected: {quinazoline_formation_detected}")

    return quinazoline_formation_detected
