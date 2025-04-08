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
    Detects if the synthesis follows a linear strategy (each step builds upon a single precursor)
    rather than a convergent strategy (combining multiple complex fragments).

    A linear strategy typically has:
    1. Most reactions with only one complex reactant
    2. Any convergent steps are early in the synthesis (high depth)
    3. Common coupling reactions may have two complex reactants but are still part of linear strategies

    Returns True if the synthesis is predominantly linear, False if convergent.
    """
    is_linear = True
    convergent_steps = []

    # Common reagents that shouldn't count toward convergence
    common_reagents = [
        "Boronic acid",
        "Boronic ester",
        "Triflate",
        "Magnesium halide",
        "Zinc halide",
        "Tin",
        "Alkyl lithium",
        "Aryl lithium",
        "Tosylate",
        "Mesylate",
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Aromatic halide",
        "Acyl halide",
        "Sulfonyl halide",
    ]

    # Coupling reactions that are common in linear synthesis
    linear_friendly_reactions = [
        "Suzuki coupling",
        "Negishi coupling",
        "Stille reaction",
        "Sonogashira",
        "Buchwald-Hartwig",
        "Heck",
        "Chan-Lam",
        "Hiyama-Denmark Coupling",
        "Kumada cross-coupling",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, convergent_steps

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a common coupling reaction
            is_coupling_reaction = False
            for reaction_type in linear_friendly_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Detected {reaction_type} at depth {depth}")
                    is_coupling_reaction = True
                    break

            # Count non-trivial reactants (more than 7 heavy atoms)
            complex_reactants = []
            for r in reactants:
                if not r:
                    continue

                mol = Chem.MolFromSmiles(r)
                if not mol:
                    continue

                # Skip if this is a common reagent
                is_common_reagent = False
                for reagent_type in common_reagents:
                    try:
                        if checker.check_fg(reagent_type, r):
                            print(f"Detected common reagent ({reagent_type}) at depth {depth}")
                            is_common_reagent = True
                            break
                    except Exception as e:
                        print(f"Error checking for {reagent_type}: {e}")

                if is_common_reagent:
                    continue

                if mol.GetNumHeavyAtoms() > 7:  # Better threshold for "complex"
                    complex_reactants.append(r)

            # If more than one complex reactant and not a coupling reaction, it's likely a convergent step
            if len(complex_reactants) > 1 and not is_coupling_reaction:
                convergent_steps.append((depth, len(complex_reactants)))
                print(
                    f"Convergent synthesis step detected at depth {depth} with {len(complex_reactants)} complex reactants"
                )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Analyze convergent steps - weigh late-stage convergence more heavily
    if convergent_steps:
        # Sort by depth (ascending - lower depth means later stage)
        convergent_steps.sort()

        # If there's late-stage convergence (depth < 4), mark as convergent
        if convergent_steps and convergent_steps[0][0] < 4:
            is_linear = False
            print(f"Late-stage convergence detected at depth {convergent_steps[0][0]}")
        # If there are multiple convergent steps, likely convergent strategy
        elif (
            len(convergent_steps) > 2
        ):  # Allow for a couple of convergent steps in linear synthesis
            is_linear = False
            print(f"Multiple convergent steps detected: {len(convergent_steps)}")

    if is_linear:
        print("Linear synthesis strategy detected")
    else:
        print("Convergent synthesis strategy detected")

    return is_linear
