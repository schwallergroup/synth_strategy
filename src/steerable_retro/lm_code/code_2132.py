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
    Detects if a cyclopropane motif is present in the final product and maintained
    throughout the synthesis.
    """
    # Check if final product has cyclopropane
    if route["type"] == "mol":
        final_product_has_cyclopropane = checker.check_ring("cyclopropane", route["smiles"])
        print(f"Final product has cyclopropane: {final_product_has_cyclopropane}")

        if not final_product_has_cyclopropane:
            return False
    else:
        print("Error: Root node is not a molecule")
        return False

    # Track if cyclopropane is maintained throughout synthesis
    cyclopropane_broken = False

    def dfs_traverse(node, depth=0):
        nonlocal cyclopropane_broken

        if cyclopropane_broken:
            return  # Stop traversal if we already found a break

        if node["type"] == "mol":
            # Skip checking starting materials
            if node.get("in_stock", False):
                return

            # For non-starting materials, check if they have cyclopropane
            has_cyclopropane = checker.check_ring("cyclopropane", node["smiles"])
            print(f"Molecule at depth {depth} has cyclopropane: {has_cyclopropane}")

        elif node["type"] == "reaction":
            # Check if reaction maintains cyclopropane
            try:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]
                reactants = rsmi.split(">")[0].split(".")

                product_has_cyclopropane = checker.check_ring("cyclopropane", product)

                # Check which reactants have cyclopropane
                reactants_with_cyclopropane = [
                    r for r in reactants if checker.check_ring("cyclopropane", r)
                ]
                reactant_has_cyclopropane = len(reactants_with_cyclopropane) > 0

                print(
                    f"Reaction at depth {depth}: product has cyclopropane: {product_has_cyclopropane}, reactants with cyclopropane: {len(reactants_with_cyclopropane)}"
                )

                # In retrosynthesis:
                # If product has cyclopropane but no reactant does, cyclopropane is formed in this step (forward)
                # This is fine - we're looking for maintenance, not formation
                if product_has_cyclopropane and not reactant_has_cyclopropane:
                    print(f"Cyclopropane formed in reaction at depth {depth} (forward direction)")

                # If product doesn't have cyclopropane but a reactant does, it's broken in this step (forward)
                elif not product_has_cyclopropane and reactant_has_cyclopropane:
                    cyclopropane_broken = True
                    print(f"Cyclopropane broken in reaction at depth {depth} (forward direction)")

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # If final product has cyclopropane and it's never broken, return True
    result = final_product_has_cyclopropane and not cyclopropane_broken
    print(f"Cyclopropane maintained: {result}")
    return result
