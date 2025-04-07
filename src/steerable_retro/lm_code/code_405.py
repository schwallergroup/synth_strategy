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
    Detects a strategy involving early-stage C-C bond formation via cross-coupling reactions
    such as Sonogashira coupling, Suzuki coupling, Negishi coupling, etc.
    """
    found_cc_coupling = False

    # List of C-C coupling reactions to check
    cc_coupling_reactions = [
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_aryl halide",
        "Sonogashira acetylene_aryl OTf",
        "Sonogashira alkyne_aryl OTf",
        "Sonogashira acetylene_alkenyl halide",
        "Sonogashira alkyne_alkenyl halide",
        "Sonogashira acetylene_alkenyl OTf",
        "Sonogashira alkyne_alkenyl OTf",
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Negishi coupling",
        "Stille reaction_aryl",
        "Stille reaction_vinyl",
        "Heck terminal vinyl",
        "Kumada cross-coupling",
        "Hiyama-Denmark Coupling",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_cc_coupling

        if node["type"] == "reaction" and depth >= 3:  # Focus on early-stage reactions
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Check if this is a C-C coupling reaction
                for reaction_type in cc_coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found C-C coupling reaction: {reaction_type} at depth {depth}")
                        print(f"Reaction SMILES: {rsmi}")
                        found_cc_coupling = True
                        break

                # If no specific reaction type matched, check for general C-C coupling patterns
                if not found_cc_coupling:
                    # Check for presence of common functional groups in C-C couplings
                    reactants_str = rsmi.split(">")[0]
                    product_str = rsmi.split(">")[-1]

                    # Check for aryl halides, terminal alkynes, boronic acids/esters
                    has_aryl_halide = checker.check_fg("Aromatic halide", reactants_str)
                    has_terminal_alkyne = checker.check_fg("Alkyne", reactants_str)
                    has_boronic = checker.check_fg(
                        "Boronic acid", reactants_str
                    ) or checker.check_fg("Boronic ester", reactants_str)

                    # Check if product has a new C-C bond (simplified check)
                    if has_aryl_halide and (has_terminal_alkyne or has_boronic):
                        product_mol = Chem.MolFromSmiles(product_str)
                        if product_mol:
                            # This is a simplified check - in a real scenario, we would
                            # need to track atom mappings to verify the new C-C bond
                            print(f"Found potential C-C coupling reaction at depth {depth}")
                            print(f"Reaction SMILES: {rsmi}")
                            found_cc_coupling = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_cc_coupling
