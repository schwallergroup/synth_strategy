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
    Detects a linear synthesis strategy with multiple fragment coupling steps.
    """
    # Track fragment coupling reactions
    fragment_coupling_reactions = 0
    linear_synthesis = True

    # List of common coupling reaction types
    coupling_reaction_types = [
        "Suzuki",
        "Negishi",
        "Stille",
        "Heck",
        "Sonogashira",
        "Buchwald-Hartwig",
        "Ullmann-Goldberg",
        "N-arylation",
    ]

    # Track the main synthesis path
    main_path_nodes = {}

    def dfs_traverse(node, depth=0):
        nonlocal fragment_coupling_reactions, linear_synthesis, main_path_nodes

        # For molecule nodes that aren't starting materials, track the main path
        if node["type"] == "mol" and not node.get("in_stock", False):
            if depth in main_path_nodes:
                # If we've already seen a different molecule at this depth, we have a branching path
                if main_path_nodes[depth]["smiles"] != node["smiles"]:
                    linear_synthesis = False
                    print(
                        f"Branching path detected at depth {depth}: {node['smiles']} vs {main_path_nodes[depth]['smiles']}"
                    )
            else:
                main_path_nodes[depth] = node

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Check if this is a fragment coupling (more than one significant reactant)
            significant_reactants = [
                r for r in reactants_smiles if len(r) > 10
            ]  # Filter out small molecules/reagents

            if len(significant_reactants) > 1:
                # Check if this is a known coupling reaction
                is_coupling = False
                for rxn_type in coupling_reaction_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_coupling = True
                        print(f"Detected {rxn_type} coupling reaction")
                        break

                # If it's a significant reactant coupling or a known coupling reaction
                if is_coupling or len([r for r in significant_reactants if len(r) > 20]) > 1:
                    fragment_coupling_reactions += 1
                    print(f"Fragment coupling detected: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if it's a linear synthesis with at least 2 fragment couplings
    strategy_present = linear_synthesis and fragment_coupling_reactions >= 2

    print(f"Linear synthesis with fragment coupling strategy:")
    print(f"  Linear synthesis: {linear_synthesis}")
    print(f"  Fragment coupling reactions: {fragment_coupling_reactions}")
    print(f"  Strategy present: {strategy_present}")

    return strategy_present
