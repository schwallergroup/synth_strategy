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
    This function detects preservation of a nitrile group throughout the synthesis.
    """
    # Track molecules with nitrile groups
    nitrile_molecules = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check if molecule contains nitrile group
            if checker.check_fg("Nitrile", mol_smiles):
                print(f"Detected nitrile in molecule at depth {depth}: {mol_smiles}")
                # Store molecule with its depth
                nitrile_molecules.append((mol_smiles, depth, node))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort molecules by depth (ascending - final product first)
    nitrile_molecules.sort(key=lambda x: x[1])

    # Check if we have nitrile in final product and at least one intermediate
    if not nitrile_molecules:
        print("No nitrile groups found in synthesis")
        return False

    # Check if final product has nitrile (should be at lowest depth)
    final_product = nitrile_molecules[0]
    if final_product[1] != 0:
        print("Final product does not contain nitrile")
        return False

    # Check if we have at least one intermediate with nitrile
    if len(nitrile_molecules) < 2:
        print("Nitrile only found in final product, not in intermediates")
        return False

    # Check if nitrile is preserved through reactions
    # This requires analyzing reaction nodes between nitrile-containing molecules

    # Get the final product node
    final_product_node = final_product[2]

    # Check if there's a path of reactions preserving the nitrile
    nitrile_preserved = False

    # For a simple check, verify if at least one intermediate and the final product have nitrile
    intermediates_with_nitrile = [m for m in nitrile_molecules if m[1] > 0]
    if intermediates_with_nitrile:
        nitrile_preserved = True
        print("Nitrile preservation strategy detected")

    return nitrile_preserved
