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
    Detects if the synthesis uses a convergent approach where two nitrile-containing
    fragments are combined in a mid-stage reaction.
    """
    has_convergent_nitrile_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_convergent_nitrile_coupling

        if node["type"] == "reaction" and depth > 0:  # Any non-terminal reaction
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if we have multiple reactants with nitriles
                nitrile_reactants = []
                for reactant in reactants_smiles:
                    if checker.check_fg("Nitrile", reactant):
                        nitrile_reactants.append(reactant)

                # Verify this is a convergent synthesis (multiple nitrile fragments combining)
                if len(nitrile_reactants) >= 2:
                    # This is a convergent nitrile synthesis - two nitrile fragments are combined
                    has_convergent_nitrile_coupling = True
                    print(f"Convergent coupling of nitrile fragments at depth {depth}")
                    print(f"Nitrile-containing reactants: {nitrile_reactants}")
                    print(f"Product: {product_smiles}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_convergent_nitrile_coupling
