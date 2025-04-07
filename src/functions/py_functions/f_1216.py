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
    """Check for nitro reduction to amine in the synthesis route"""
    found = False

    def dfs(node, depth=0):
        nonlocal found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for nitro reduction reaction
            if checker.check_reaction(
                "Reduction of nitro groups to amines", rxn_smiles
            ):
                found = True
                print(f"Found nitro reduction at depth {depth}: {rxn_smiles}")

                # Verify that nitro group is present in reactants and amine in product
                try:
                    reactants = rxn_smiles.split(">")[0].split(".")
                    product = rxn_smiles.split(">")[-1]

                    has_nitro_reactant = any(
                        checker.check_fg("Nitro group", r) for r in reactants
                    )
                    has_amine_product = (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                        or checker.check_fg("Aniline", product)
                    )

                    found = has_nitro_reactant and has_amine_product
                except Exception as e:
                    print(f"Error checking nitro reduction: {e}")

            # Additional check for nitro to amine conversion even if not explicitly labeled
            if not found:
                try:
                    reactants = rxn_smiles.split(">")[0].split(".")
                    product = rxn_smiles.split(">")[-1]

                    has_nitro_reactant = any(
                        checker.check_fg("Nitro group", r) for r in reactants
                    )
                    has_amine_product = (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                        or checker.check_fg("Aniline", product)
                    )

                    if has_nitro_reactant and has_amine_product:
                        found = True
                        print(
                            f"Found implicit nitro reduction at depth {depth}: {rxn_smiles}"
                        )
                except Exception as e:
                    print(f"Error checking implicit nitro reduction: {e}")

        # Recursively check children
        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)
    return found
