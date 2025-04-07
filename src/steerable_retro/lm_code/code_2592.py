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
    This function detects a synthetic strategy involving late-stage reductive amination
    (aldehyde/ketone + amine → amine bond formation in the final steps).
    """
    found_reductive_amination = False

    def dfs_traverse(node, depth=0):
        nonlocal found_reductive_amination

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a reductive amination reaction
            is_reductive_amination = (
                checker.check_reaction("Reductive amination with aldehyde", rsmi)
                or checker.check_reaction("Reductive amination with ketone", rsmi)
                or checker.check_reaction("Reductive amination with alcohol", rsmi)
            )

            # If direct reaction check fails, try checking by functional groups
            if not is_reductive_amination:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                has_carbonyl = False
                has_amine = False

                for reactant in reactants:
                    if not reactant:
                        continue

                    # Check for aldehyde or ketone in reactants
                    if checker.check_fg("Aldehyde", reactant) or checker.check_fg(
                        "Ketone", reactant
                    ):
                        has_carbonyl = True
                        print(f"Found carbonyl group in reactant: {reactant}")

                    # Check for primary or secondary amine in reactants
                    if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                        "Secondary amine", reactant
                    ):
                        has_amine = True
                        print(f"Found amine group in reactant: {reactant}")

                # Check if product has a new C-N bond that wasn't in reactants
                if has_carbonyl and has_amine:
                    product_mol = Chem.MolFromSmiles(product)
                    if (
                        product_mol
                        and checker.check_fg("Tertiary amine", product)
                        or checker.check_fg("Secondary amine", product)
                    ):
                        print(f"Found potential reductive amination at depth {depth}: {rsmi}")
                        is_reductive_amination = True

            # If this is a late-stage reaction (depth <= 1) and is reductive amination
            if depth <= 1 and is_reductive_amination:
                found_reductive_amination = True
                print(f"Detected late-stage reductive amination at depth {depth}: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Late-stage reductive amination strategy detected: {found_reductive_amination}")
    return found_reductive_amination
