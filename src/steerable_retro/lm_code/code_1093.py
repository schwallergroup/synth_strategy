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
    This function detects if the synthetic route maintains a cyano group
    throughout the synthesis.
    """
    # Track the final product (target molecule)
    final_product_smiles = route["smiles"]
    print(f"Final product: {final_product_smiles}")

    # Check if final product has cyano group
    if not checker.check_fg("Nitrile", final_product_smiles):
        print(f"Final product does not have a cyano group: {final_product_smiles}")
        return False

    # Track molecules in the synthetic pathway that should maintain the cyano group
    cyano_maintained = True

    def dfs_traverse(node, depth=0, in_cyano_path=False):
        nonlocal cyano_maintained

        # For molecule nodes
        if node["type"] == "mol" and node.get("smiles"):
            mol_smiles = node["smiles"]
            has_cyano = checker.check_fg("Nitrile", mol_smiles)

            # If this is the first node (final product), mark it as in the cyano path
            if depth == 0:
                in_cyano_path = has_cyano
                print(f"Depth {depth}, Final product has cyano: {has_cyano}")

            # If we're in a path that should maintain cyano but this molecule doesn't have it
            if in_cyano_path and not has_cyano and not node.get("in_stock", False):
                print(f"Molecule at depth {depth} without cyano group found: {mol_smiles}")
                cyano_maintained = False

        # For reaction nodes, check if the cyano group is maintained through the reaction
        elif (
            node["type"] == "reaction"
            and in_cyano_path
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]
            reactants = rsmi.split(">")[0].split(".")

            # Check if product has cyano
            product_has_cyano = checker.check_fg("Nitrile", product)

            # Check if at least one reactant has cyano
            reactant_has_cyano = any(checker.check_fg("Nitrile", r) for r in reactants)

            # If cyano is in product but not in any reactant, it's being introduced
            if product_has_cyano and not reactant_has_cyano:
                print(f"Cyano group introduced in reaction: {rsmi}")
                # This is fine - cyano is being introduced

            # If cyano is in reactant but not in product, it's being removed
            elif reactant_has_cyano and not product_has_cyano:
                print(f"Cyano group removed in reaction: {rsmi}")
                cyano_maintained = False

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, in_cyano_path)

    # Start traversal from root
    dfs_traverse(route)
    return cyano_maintained
