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
    Checks if the route contains late-stage sulfonamide formation.
    Late-stage means the sulfonamide is formed in the final steps (low depth).
    """
    # Track sulfonamide formation reactions and their depths
    sulfonamide_reactions = []
    max_depth = 0

    def dfs(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        # Check for sulfonamide formation in reaction nodes
        if (
            node.get("type") == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check if the product contains a sulfonamide group
            if checker.check_fg("Sulfonamide", product):
                # Check if reactants don't have sulfonamide (indicating formation)
                reactants = reactants_part.split(".")
                has_sulfonamide_in_reactants = False

                for reactant in reactants:
                    print(f"Checking reactant for sulfonamide: {reactant}")
                    if checker.check_fg("Sulfonamide", reactant):
                        has_sulfonamide_in_reactants = True
                        print(f"Found sulfonamide in reactant: {reactant}")
                        break

                print(
                    f"Product has sulfonamide: {checker.check_fg('Sulfonamide', product)}"
                )
                print(f"Any reactant has sulfonamide: {has_sulfonamide_in_reactants}")

                # If product has sulfonamide but reactants don't, it's a formation reaction
                if not has_sulfonamide_in_reactants:
                    print(f"Found sulfonamide formation at depth {depth}: {rsmi}")
                    sulfonamide_reactions.append((depth, rsmi))

                    # Also check for specific sulfonamide formation reactions
                    is_schotten_baumann = checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                    ) or checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                    )
                    print(
                        f"Is Schotten-Baumann sulfonamide synthesis: {is_schotten_baumann}"
                    )

        # Recursively check children
        for child in node.get("children", []):
            dfs(child, depth + 1)

    # Start DFS from the root
    dfs(route)

    print(f"Max depth in synthesis: {max_depth}")
    print(f"Sulfonamide reactions found: {len(sulfonamide_reactions)}")

    # If we found sulfonamide formations, check if any are late-stage
    if sulfonamide_reactions:
        # Sort by depth (ascending)
        sulfonamide_reactions.sort(key=lambda x: x[0])
        print(
            f"Sulfonamide reaction depths: {[depth for depth, _ in sulfonamide_reactions]}"
        )

        # Consider it late-stage if it's in the first third of the synthesis
        late_stage_threshold = max(2, max_depth // 3)
        print(f"Late-stage threshold: {late_stage_threshold}")

        return sulfonamide_reactions[0][0] <= late_stage_threshold

    return False
