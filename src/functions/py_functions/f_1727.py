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
    This function detects a strategy where a multi-component reaction (MCR)
    is used in the final step of the synthesis.
    """
    # Initialize tracking
    has_late_stage_mcr = False

    # Define what constitutes "late stage" - typically first 2 levels of the tree
    LATE_STAGE_DEPTH_THRESHOLD = 2

    # List of known MCR reaction types
    mcr_reaction_types = ["Ugi reaction", "A3 coupling", "A3 coupling to imidazoles"]

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_mcr

        # Check if this is a reaction node at a late stage (low depth)
        if node["type"] == "reaction" and depth <= LATE_STAGE_DEPTH_THRESHOLD:
            # Extract reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            # Extract reactants
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Check if it has 3+ components (potential MCR)
            if len(reactants_smiles) >= 3:
                print(
                    f"Found potential MCR with {len(reactants_smiles)} components at depth {depth}"
                )

                # Check if it matches known MCR reaction types
                for mcr_type in mcr_reaction_types:
                    if checker.check_reaction(mcr_type, rsmi):
                        print(f"Confirmed as {mcr_type}")
                        has_late_stage_mcr = True
                        return

                # If no specific MCR type matched but it has 3+ components,
                # we'll still consider it an MCR for this purpose
                has_late_stage_mcr = True

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root (final product)
    dfs_traverse(route)

    return has_late_stage_mcr
