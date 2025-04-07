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
    This function detects a linear synthesis strategy with late-stage fragment coupling.
    """
    # Track if we found a linear synthesis with late coupling
    found_linear_with_late_coupling = False

    # Track the depth of the first coupling step
    first_coupling_depth = None
    max_depth = 0

    # Track branching at each depth to check linearity
    reactions_at_depth = {}

    # List of coupling reaction types to check
    coupling_reaction_types = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Suzuki coupling with boronic acids OTf",
        "Suzuki coupling with boronic esters OTf",
        "Negishi coupling",
        "Stille reaction_vinyl",
        "Stille reaction_aryl",
        "Stille reaction_benzyl",
        "Stille reaction_allyl",
        "Stille reaction_vinyl OTf",
        "Stille reaction_aryl OTf",
        "Stille reaction_benzyl OTf",
        "Stille reaction_allyl OTf",
        "Stille reaction_other",
        "Stille reaction_other OTf",
        "Heck terminal vinyl",
        "Heck_terminal_vinyl",
        "Heck_non-terminal_vinyl",
        "Hiyama-Denmark Coupling",
        "Kumada cross-coupling",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Ullmann-Goldberg Substitution amine",
        "Ullmann-Goldberg Substitution thiol",
        "Ullmann-Goldberg Substitution aryl alcohol",
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_aryl halide",
        "Sonogashira acetylene_aryl OTf",
        "Sonogashira alkyne_aryl OTf",
        "Sonogashira acetylene_alkenyl halide",
        "Sonogashira alkyne_alkenyl halide",
        "Sonogashira acetylene_alkenyl OTf",
        "Sonogashira alkyne_alkenyl OTf",
        "Sonogashira acetylene_acyl halide",
        "Sonogashira alkyne_acyl halide",
    ]

    def is_coupling_reaction(rsmi):
        """Check if a reaction is a coupling reaction"""
        # First check if it has multiple reactants
        reactants = rsmi.split(">")[0].split(".")
        if len(reactants) <= 1:
            return False

        # Then check if it matches any known coupling reaction types
        for rxn_type in coupling_reaction_types:
            if checker.check_reaction(rxn_type, rsmi):
                print(f"Found coupling reaction: {rxn_type}")
                return True

        # If no specific coupling reaction type matched, check for multiple reactants
        # and at least one aromatic or alkene/alkyne component
        has_aromatic = False
        has_alkene_alkyne = False

        for reactant in reactants:
            mol = Chem.MolFromSmiles(reactant)
            if mol:
                # Check for aromatic atoms
                for atom in mol.GetAtoms():
                    if atom.GetIsAromatic():
                        has_aromatic = True
                        break

                # Check for alkene/alkyne
                if (
                    checker.check_fg("Alkyne", reactant)
                    or checker.check_fg("Vinyl", reactant)
                    or checker.check_fg("Allyl", reactant)
                ):
                    has_alkene_alkyne = True

        # If we have multiple reactants and at least one has aromatic or alkene/alkyne character,
        # it might be a coupling reaction
        return has_aromatic or has_alkene_alkyne

    def dfs_traverse(node, depth=0):
        nonlocal found_linear_with_late_coupling, first_coupling_depth, max_depth

        # Update max depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Track number of reactions at this depth
            if depth not in reactions_at_depth:
                reactions_at_depth[depth] = 0
            reactions_at_depth[depth] += 1

            # Check if this is a coupling reaction
            if is_coupling_reaction(rsmi):
                if first_coupling_depth is None or depth < first_coupling_depth:
                    first_coupling_depth = depth
                    print(f"Found coupling reaction at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if synthesis is linear before the coupling
    is_linear = True
    if first_coupling_depth is not None:
        for d in range(first_coupling_depth + 1, max_depth + 1):
            if d in reactions_at_depth and reactions_at_depth[d] > 1:
                is_linear = False
                break

    # If the first coupling is at depth 0 or 1 (late stage), we have a reasonable synthesis depth,
    # and the synthesis is linear after the coupling
    if (
        first_coupling_depth is not None
        and first_coupling_depth <= 1
        and max_depth >= 3
        and is_linear
    ):
        found_linear_with_late_coupling = True
        print(f"Found linear synthesis with late-stage coupling at depth {first_coupling_depth}")
        print(f"Max depth: {max_depth}")
        print(f"Reactions at each depth: {reactions_at_depth}")

    return found_linear_with_late_coupling
