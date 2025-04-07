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
    Detects if the synthesis involves protection of a ketone as a cyclic ketal.

    In retrosynthetic analysis, we're looking for ketal deprotection to ketone,
    which corresponds to ketone protection in the forward direction.
    """
    ketone_protection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal ketone_protection_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if this is a ketal deprotection reaction (retrosynthetically)
                    if checker.check_reaction("Ketal hydrolysis to ketone", rsmi):
                        print(f"Ketal hydrolysis to ketone reaction detected at depth {depth}")
                        ketone_protection_found = True
                    elif checker.check_reaction("Acetal hydrolysis to ketone", rsmi):
                        print(f"Acetal hydrolysis to ketone reaction detected at depth {depth}")
                        ketone_protection_found = True
                    else:
                        # Manual check for ketone protection/deprotection
                        product_has_ketone = checker.check_fg("Ketone", product)

                        # Check if any reactant has a cyclic ketal structure
                        reactant_has_ketal = False
                        for reactant in reactants:
                            # Check for dioxolane or dioxane rings which are common ketal protecting groups
                            if checker.check_ring("dioxolane", reactant) or checker.check_ring(
                                "dioxane", reactant
                            ):
                                reactant_has_ketal = True
                                break

                        if product_has_ketone and reactant_has_ketal:
                            print(
                                f"Ketone deprotection detected at depth {depth} (product has ketone, reactant has ketal)"
                            )
                            ketone_protection_found = True

                        # Also check the forward direction (for completeness)
                        reactant_has_ketone = any(checker.check_fg("Ketone", r) for r in reactants)
                        product_has_ketal = checker.check_ring(
                            "dioxolane", product
                        ) or checker.check_ring("dioxane", product)

                        if (
                            reactant_has_ketone
                            and product_has_ketal
                            and checker.check_reaction("Aldehyde or ketone acetalization", rsmi)
                        ):
                            print(
                                f"Ketone protection detected at depth {depth} (reactant has ketone, product has ketal)"
                            )
                            ketone_protection_found = True

                except Exception as e:
                    print(f"Error processing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return ketone_protection_found
