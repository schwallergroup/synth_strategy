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
    This function detects a synthetic strategy involving ring opening,
    specifically ketal/acetal ring opening and other ring opening reactions.
    """
    has_ring_opening = False

    def dfs_traverse(node):
        nonlocal has_ring_opening

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for specific ring opening reactions
            if (
                checker.check_reaction("Acetal hydrolysis to aldehyde", rsmi)
                or checker.check_reaction("Ketal hydrolysis to ketone", rsmi)
                or checker.check_reaction("Acetal hydrolysis to diol", rsmi)
            ):
                print(f"Found ketal/acetal ring opening reaction: {rsmi}")
                has_ring_opening = True
                return

            # Check for epoxide ring opening
            if checker.check_reaction("Ring opening of epoxide with amine", rsmi):
                print(f"Found epoxide ring opening reaction: {rsmi}")
                has_ring_opening = True
                return

            # Check for other ring opening patterns
            for reactant in reactants_smiles:
                # Check for cyclic structures in reactants
                has_cyclic_reactant = any(
                    checker.check_ring(ring_name, reactant)
                    for ring_name in [
                        "dioxane",
                        "dioxolane",
                        "oxirane",
                        "oxetane",
                        "tetrahydrofuran",
                        "tetrahydropyran",
                        "aziridine",
                        "azetidine",
                        "pyrrolidine",
                        "piperidine",
                    ]
                )

                if has_cyclic_reactant:
                    # Try to confirm by counting rings
                    try:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        product_mol = Chem.MolFromSmiles(product_smiles)

                        if reactant_mol and product_mol:
                            reactant_rings = Chem.GetSSSR(reactant_mol)
                            product_rings = Chem.GetSSSR(product_mol)

                            if len(reactant_rings) > len(product_rings):
                                print(f"Found ring opening based on ring count: {rsmi}")
                                print(
                                    f"Reactant rings: {len(reactant_rings)}, Product rings: {len(product_rings)}"
                                )
                                has_ring_opening = True
                                return
                    except Exception as e:
                        print(f"Error in ring counting: {e}")

                    # Check for functional group changes consistent with ring opening
                    if (
                        checker.check_ring("dioxane", reactant)
                        or checker.check_ring("dioxolane", reactant)
                    ) and (
                        checker.check_fg("Ketone", product_smiles)
                        or checker.check_fg("Aldehyde", product_smiles)
                    ):
                        print(f"Found ketal/acetal ring opening based on functional groups: {rsmi}")
                        has_ring_opening = True
                        return

                    if checker.check_ring("oxirane", reactant) and "OH" in product_smiles:
                        print(f"Found epoxide ring opening based on functional groups: {rsmi}")
                        has_ring_opening = True
                        return

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    if has_ring_opening:
        print("Detected ring opening strategy")
    return has_ring_opening
