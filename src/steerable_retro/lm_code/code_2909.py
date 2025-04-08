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
    This function detects a strategy where a nitro group is preserved throughout
    the entire synthesis.
    """
    # First check if the target molecule has a nitro group
    if route["type"] == "mol" and "smiles" in route:
        target_mol_smiles = route["smiles"]
        if not checker.check_fg("Nitro group", target_mol_smiles):
            print(f"Target molecule does not have a nitro group: {target_mol_smiles}")
            return False
        print(f"Target molecule has a nitro group: {target_mol_smiles}")
    else:
        print("Invalid route: Root node is not a molecule")
        return False

    # Track molecules in the synthetic pathway
    synthetic_pathway_molecules = []

    def dfs_traverse(node, depth=0):
        # For molecule nodes that aren't starting materials
        if node["type"] == "mol" and "smiles" in node and not node.get("in_stock", False):
            synthetic_pathway_molecules.append((node["smiles"], depth))

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort molecules by depth (from target to starting materials)
    synthetic_pathway_molecules.sort(key=lambda x: x[1])

    # Check if all intermediate molecules have a nitro group
    nitro_count = 0
    total_count = len(synthetic_pathway_molecules)

    for smiles, depth in synthetic_pathway_molecules:
        if checker.check_fg("Nitro group", smiles):
            nitro_count += 1
            print(f"Found nitro group in intermediate at depth {depth}: {smiles}")
        else:
            print(f"Missing nitro group in intermediate at depth {depth}: {smiles}")

    # Return True if all molecules in the synthetic pathway have a nitro group
    # and there's at least one intermediate molecule
    result = total_count > 0 and nitro_count == total_count
    print(
        f"Nitro group preservation detected: {result} (nitro count: {nitro_count}, total intermediates: {total_count})"
    )
    return result
