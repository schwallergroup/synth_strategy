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
    This function detects if a cyano group is maintained throughout the synthesis route.

    Returns True if a cyano group is present in the target molecule and persists
    throughout the entire synthesis route (i.e., is present in all molecules and
    preserved across all reactions).
    """
    # Check if target molecule has a cyano group
    if not route or route["type"] != "mol" or not route["smiles"]:
        print("Invalid route or no target molecule")
        return False

    target_has_cyano = checker.check_fg("Nitrile", route["smiles"])
    if not target_has_cyano:
        print(f"Target molecule has no cyano group: {route['smiles']}")
        return False

    # Get cyano atom indices in target molecule
    target_cyano_indices = checker.get_fg_atom_indices("Nitrile", route["smiles"])
    if not target_cyano_indices:
        print(f"Could not identify cyano group atoms in target: {route['smiles']}")
        return False

    print(f"Target molecule has cyano group(s): {target_cyano_indices}")

    # Track if cyano group persists through non-starting material molecules
    cyano_persists = True

    # Keep track of atom mappings for cyano groups
    cyano_atom_mappings = {}

    def dfs_traverse(node, depth=0):
        nonlocal cyano_persists

        if not cyano_persists:  # Early termination if we already found an issue
            return

        if node["type"] == "mol" and node["smiles"]:
            # Skip check for starting materials
            if node.get("in_stock", False):
                print(f"Depth {depth}: Starting material (skipping cyano check): {node['smiles']}")
                return

            # Check if molecule has a cyano group
            has_cyano = checker.check_fg("Nitrile", node["smiles"])
            if not has_cyano:
                cyano_persists = False
                print(f"Depth {depth}: Molecule without cyano group found: {node['smiles']}")
                return
            print(f"Depth {depth}: Molecule has cyano group: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Depth {depth}: Processing reaction: {rsmi}")

            # Check if reaction preserves cyano group
            # For reactions, we need to verify the cyano group is preserved
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has cyano group
                if not checker.check_fg("Nitrile", product):
                    cyano_persists = False
                    print(f"Depth {depth}: Reaction product doesn't have cyano group: {product}")
                    return

                # Check if at least one reactant has cyano group
                reactant_has_cyano = any(checker.check_fg("Nitrile", r) for r in reactants)
                if not reactant_has_cyano:
                    cyano_persists = False
                    print(f"Depth {depth}: No reactant has cyano group, but product does: {rsmi}")
                    return
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    return cyano_persists
