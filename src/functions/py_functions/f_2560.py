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
    This function detects a linear synthesis strategy with a late-stage coupling
    to introduce structural complexity.
    """
    late_coupling = False
    linear_steps = 0
    coupling_reaction_types = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "Heck terminal vinyl",
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_aryl halide",
        "Stille reaction_aryl",
        "Negishi coupling",
        "Ullmann condensation",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal late_coupling, linear_steps

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")

            # Check if it's a coupling reaction (at late stage)
            if depth <= 2 and len(reactants) >= 2:
                # Check for complex reactants
                complex_reactants = sum(
                    1
                    for r in reactants
                    if Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).GetNumAtoms() > 8
                )

                # Check if it's a known coupling reaction type
                is_coupling_reaction = False
                for rxn_type in coupling_reaction_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_coupling_reaction = True
                        print(
                            f"Found late-stage coupling reaction: {rxn_type} at depth {depth}"
                        )
                        break

                if (complex_reactants >= 2 and is_coupling_reaction) or (
                    complex_reactants >= 2 and len(reactants) == 2
                ):
                    late_coupling = True
                    print(
                        f"Confirmed late-stage coupling at depth {depth} with {complex_reactants} complex reactants"
                    )

            # Count linear steps (single complex reactant)
            if len(reactants) >= 1:
                complex_reactants = sum(
                    1
                    for r in reactants
                    if Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).GetNumAtoms() > 8
                )
                if complex_reactants == 1 and len(reactants) <= 2:
                    linear_steps += 1
                    print(f"Found linear synthesis step at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Summary: late_coupling={late_coupling}, linear_steps={linear_steps}")
    # Return True if we found a late-stage coupling and at least 3 linear steps
    return late_coupling and linear_steps >= 3
