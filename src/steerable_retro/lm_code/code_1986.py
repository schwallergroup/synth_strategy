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
    Detects if the synthesis route involves borylation of an aryl compound for cross-coupling.
    """
    # Track borylation reactions and their products
    borylation_reactions = []

    def dfs_traverse(node, depth=0, path=None):
        if path is None:
            path = []

        current_path = path + [node]

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for borylation reactions
            is_borylation = False

            # Check for specific borylation reaction types
            if (
                checker.check_reaction("Preparation of boronic acids", rsmi)
                or checker.check_reaction(
                    "Preparation of boronic acids without boronic ether", rsmi
                )
                or checker.check_reaction(
                    "Preparation of boronic acids from trifluoroborates", rsmi
                )
                or checker.check_reaction("Preparation of boronic ethers", rsmi)
            ):
                is_borylation = True
                print(f"Detected borylation reaction at depth {depth}: {rsmi}")

            # Alternative check: product has boronic acid/ester but reactants don't
            if not is_borylation:
                has_boronic_in_product = checker.check_fg(
                    "Boronic acid", product
                ) or checker.check_fg("Boronic ester", product)

                has_boronic_in_reactants = False
                has_aromatic_halide_in_reactants = False

                for reactant in reactants:
                    if checker.check_fg("Boronic acid", reactant) or checker.check_fg(
                        "Boronic ester", reactant
                    ):
                        has_boronic_in_reactants = True
                    if checker.check_fg("Aromatic halide", reactant):
                        has_aromatic_halide_in_reactants = True

                if (
                    has_boronic_in_product
                    and not has_boronic_in_reactants
                    and has_aromatic_halide_in_reactants
                ):
                    is_borylation = True
                    print(f"Detected borylation pattern at depth {depth}: {rsmi}")

            # If borylation detected, store the reaction and product for later checking
            if is_borylation:
                borylation_reactions.append((product, depth, current_path))

            # Check if this is a coupling reaction using a previously formed boronic compound
            if (
                checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
            ):

                # Check if any reactant was produced by a previous borylation
                for reactant in reactants:
                    for borylation_product, _, _ in borylation_reactions:
                        if reactant == borylation_product:
                            print(
                                f"Found coupling reaction using borylated compound at depth {depth}"
                            )
                            return True

        # Traverse children
        for child in node.get("children", []):
            if dfs_traverse(child, depth + 1, current_path):
                return True

        return False

    # Start DFS traversal
    if dfs_traverse(route):
        return True

    # If we found borylation reactions but didn't confirm coupling use,
    # check if any borylation product appears in a later reaction
    for borylation_product, depth, path in borylation_reactions:
        # Check if this product is used in any subsequent reaction
        for i, node in enumerate(path):
            if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Skip reactions that come before this borylation in the path
                if i <= depth:
                    continue

                # Check if the borylated product is used in this reaction
                for reactant in reactants:
                    if reactant == borylation_product:
                        print(f"Borylated compound used in subsequent reaction")
                        return True

    return len(borylation_reactions) > 0
