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
    Detects if the synthesis follows a linear pattern without convergent steps.

    A linear synthesis is characterized by:
    1. Each reaction step has at most 2 reactants (one main substrate and one reagent)
    2. The synthesis path doesn't branch out (each non-leaf node has exactly one child)
    """
    is_linear = True
    max_reactants_per_step = 0

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, max_reactants_per_step

        if node["type"] == "reaction":
            # Count number of molecule children (precursors in retrosynthesis)
            mol_children = [
                child for child in node.get("children", []) if child["type"] == "mol"
            ]
            reactant_count = len(mol_children)
            max_reactants_per_step = max(max_reactants_per_step, reactant_count)

            # Check if this is a known multicomponent reaction that's still considered linear
            if reactant_count > 2 and "metadata" in node and "rsmi" in node["metadata"]:
                rxn_smiles = node["metadata"]["rsmi"]
                # Check for common multicomponent reactions that are still considered linear
                if (
                    checker.check_reaction("Ugi reaction", rxn_smiles)
                    or checker.check_reaction("A3 coupling", rxn_smiles)
                    or checker.check_reaction("A3 coupling to imidazoles", rxn_smiles)
                ):
                    print(
                        f"Multicomponent reaction detected at depth {depth}, still considered linear"
                    )
                    # Don't mark as non-linear for these specific reactions
                else:
                    # Check if most reactants are simple starting materials
                    complex_intermediates = sum(
                        1 for child in mol_children if not child.get("in_stock", False)
                    )
                    if complex_intermediates > 1:
                        is_linear = False
                        print(
                            f"Detected convergent synthesis at depth {depth}: {complex_intermediates} complex intermediates converge"
                        )
                    else:
                        print(
                            f"Multiple reactants at depth {depth}, but only one complex intermediate"
                        )
            elif reactant_count > 2:
                # If we can't check reaction type, use the default logic
                is_linear = False
                print(
                    f"Detected non-linear synthesis at depth {depth}: {reactant_count} reactants in one step"
                )

            # Check if any non-terminal molecule node has multiple reaction children
            # This would indicate branching in the synthesis path
            for child in mol_children:
                if not child.get("in_stock", False):  # Not a starting material
                    reaction_children = [
                        c for c in child.get("children", []) if c["type"] == "reaction"
                    ]
                    if len(reaction_children) > 1:
                        is_linear = False
                        print(
                            f"Detected branching synthesis at depth {depth+1}: molecule splits into {len(reaction_children)} reaction paths"
                        )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # A synthesis is linear if:
    # 1. Maximum reactants per step is ≤ 2 OR it's a recognized multicomponent reaction
    # 2. No branching detected in the synthesis path
    return is_linear
