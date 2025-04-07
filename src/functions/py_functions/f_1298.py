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
    This function detects if the synthesis involves a nitro reduction step
    (converting NO2 to NH2).
    """
    nitro_reduction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Direct check for nitro reduction reaction
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(
                        f"Nitro reduction detected at depth {depth} using reaction type check"
                    )
                    nitro_reduction_found = True
                    return

                # Check for nitro group in reactants and primary amine in product
                has_nitro_reactant = any(
                    checker.check_fg("Nitro group", r) for r in reactants
                )
                has_amine_product = checker.check_fg(
                    "Primary amine", product
                ) or checker.check_fg("Aniline", product)

                if has_nitro_reactant and has_amine_product:
                    # Additional verification that it's a reduction reaction
                    # Check if the nitro group is actually converted to an amine
                    for r in reactants:
                        if checker.check_fg("Nitro group", r):
                            # If we have a nitro group in reactant and primary amine in product,
                            # and the reaction involves reduction, it's likely a nitro reduction
                            print(
                                f"Nitro reduction detected at depth {depth} using functional group analysis"
                            )
                            nitro_reduction_found = True
                            return
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return nitro_reduction_found
