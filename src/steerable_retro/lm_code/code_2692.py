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
    Detects if the route involves formation of an indazole ring system via
    hydrazine condensation with a ketone.
    """
    indazole_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal indazole_formation

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains indazole ring
                if checker.check_ring("indazole", product):
                    print(f"Found indazole ring in product at depth {depth}")

                    # Check for hydrazine and ketone in reactants
                    has_hydrazine = False
                    has_ketone = False

                    for reactant in reactants:
                        if checker.check_fg("Hydrazine", reactant):
                            has_hydrazine = True
                            print(f"Found hydrazine in reactants at depth {depth}")
                        if checker.check_fg("Ketone", reactant):
                            has_ketone = True
                            print(f"Found ketone in reactants at depth {depth}")

                    # Check if product doesn't have hydrazine (it was consumed)
                    product_has_hydrazine = checker.check_fg("Hydrazine", product)

                    # Verify this is a condensation reaction (hydrazine + ketone → indazole)
                    if has_hydrazine and has_ketone and not product_has_hydrazine:
                        # Additional check for reaction type if possible
                        if checker.check_reaction(
                            "Hydrazone formation", rsmi
                        ) or checker.check_reaction(
                            "Addition of primary amines to ketones/thiocarbonyls", rsmi
                        ):
                            indazole_formation = True
                            print(
                                f"Confirmed indazole formation via hydrazine condensation at depth {depth}"
                            )
                        else:
                            # If specific reaction check fails, still accept if we have all the right components
                            indazole_formation = True
                            print(
                                f"Detected indazole formation via hydrazine at depth {depth} (reaction type not specifically matched)"
                            )
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return indazole_formation
