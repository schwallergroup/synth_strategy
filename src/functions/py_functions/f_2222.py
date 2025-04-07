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
    This function detects a synthesis strategy where a specific functional group
    (like a brominated aromatic ring) is preserved throughout the synthesis as a spectator.
    """
    # Track molecules with brominated aromatic rings and their positions
    molecules_with_bromine = []
    reaction_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal reaction_steps

        if node["type"] == "mol" and not node.get("in_stock", False):
            # Check if molecule has a brominated aromatic ring
            if checker.check_fg("Aromatic halide", node["smiles"]):
                # Get the molecule for further analysis
                mol = Chem.MolFromSmiles(node["smiles"])
                # Store molecule info for later comparison
                molecules_with_bromine.append(
                    {"smiles": node["smiles"], "depth": depth}
                )
                print(
                    f"Found molecule with aromatic halide at depth {depth}: {node['smiles']}"
                )

        elif node["type"] == "reaction":
            reaction_steps += 1
            # Check if the reaction involves aromatic halide coupling
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                if (
                    checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi)
                    or checker.check_reaction("Buchwald-Hartwig", rsmi)
                    or checker.check_reaction("Heck terminal vinyl", rsmi)
                ):
                    print(
                        f"Found coupling reaction that likely consumes aromatic halide: {rsmi}"
                    )
                    return  # This branch doesn't preserve the spectator group

        # Process children
        if "children" in node:
            for child in node["children"]:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have enough molecules with brominated aromatic rings
    print(f"Found {len(molecules_with_bromine)} molecules with aromatic halides")
    print(f"Reaction steps: {reaction_steps}")

    # We need at least 3 molecules (including the target) with preserved brominated aromatic rings
    # and at least 2 reaction steps (to connect 3 molecules)
    if len(molecules_with_bromine) >= 3 and reaction_steps >= 2:
        # Check if the brominated aromatic is preserved across different depths
        depths = set(mol_info["depth"] for mol_info in molecules_with_bromine)
        if (
            len(depths) >= 3
        ):  # Ensure molecules at different synthesis stages have the group
            print("Brominated aromatic ring is preserved throughout synthesis")
            return True

    print("Brominated aromatic ring is not preserved throughout synthesis")
    return False
