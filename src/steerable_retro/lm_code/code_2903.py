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
    This function detects if a nitro group is preserved throughout the synthesis.
    """
    # Track nitro groups through the synthesis
    preserved_nitro = False

    def get_nitro_atom_indices(mol_smiles):
        """Get atom indices of nitro groups in a molecule"""
        return checker.get_fg_atom_indices("Nitro group", mol_smiles)

    def dfs_traverse(node, depth=0, parent_nitro_atoms=None):
        nonlocal preserved_nitro

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if this molecule has nitro groups
            has_nitro = checker.check_fg("Nitro group", mol_smiles)

            # If we're at the target molecule (depth 0), check if it has nitro groups
            if depth == 0:
                if has_nitro:
                    # Get atom indices of nitro groups in the target molecule
                    nitro_indices = get_nitro_atom_indices(mol_smiles)
                    if nitro_indices:
                        # Continue traversal with these nitro groups to track
                        for child in node.get("children", []):
                            dfs_traverse(child, depth + 1, nitro_indices)
                else:
                    # Target molecule has no nitro groups, so no preservation to check
                    preserved_nitro = False
                    return

            # For starting materials (leaf nodes)
            elif not node.get("children", []) and node.get("in_stock", False):
                # If we've tracked a nitro group all the way to a starting material,
                # it means the nitro group is preserved throughout the synthesis
                if has_nitro and parent_nitro_atoms:
                    preserved_nitro = True
                    print(f"Nitro group preserved to starting material: {mol_smiles}")

            # For intermediate molecules
            else:
                # If we're tracking nitro groups and this molecule has them
                if has_nitro and parent_nitro_atoms:
                    # Continue tracking through children
                    nitro_indices = get_nitro_atom_indices(mol_smiles)
                    for child in node.get("children", []):
                        dfs_traverse(child, depth + 1, nitro_indices)
                else:
                    # Nitro group lost at this step
                    return

        elif node["type"] == "reaction":
            # For reaction nodes, we need to check if the nitro group is preserved
            # through the reaction by examining atom mapping
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if any reactant has a nitro group
                reactant_has_nitro = any(
                    checker.check_fg("Nitro group", r) for r in reactants_smiles
                )

                # If we're tracking nitro groups and reactants have them
                if reactant_has_nitro and parent_nitro_atoms:
                    # Continue tracking through children
                    for child in node.get("children", []):
                        dfs_traverse(child, depth + 1, parent_nitro_atoms)
            except Exception as e:
                print(f"Error processing reaction node: {e}")
                # If we can't process the reaction, assume nitro group is not preserved
                return

        # For any other node type, continue traversal
        else:
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1, parent_nitro_atoms)

    # Start traversal from the root node
    dfs_traverse(route)

    if preserved_nitro:
        print("Nitro group is preserved throughout the synthesis")
    else:
        print("Nitro group is not preserved throughout the synthesis")

    return preserved_nitro
