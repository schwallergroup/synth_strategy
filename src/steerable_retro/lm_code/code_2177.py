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
    Detects if the synthesis uses a late-stage secondary alcohol oxidation strategy.
    This checks if the final reaction (depth 0 or 1) involves oxidation of a secondary alcohol to a ketone.
    """
    # Track if we found the pattern
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern
        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and depth <= 1:  # Final reaction or one step before
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction SMILES at depth {depth}: {rsmi}")

                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for both oxidation and reduction reactions
                oxidation_reaction = checker.check_reaction(
                    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
                )
                reduction_reaction = checker.check_reaction(
                    "Reduction of ketone to secondary alcohol", rsmi
                )

                if oxidation_reaction:
                    print("Found alcohol oxidation reaction")

                    # In forward direction: alcohol → ketone
                    for reactant in reactants:
                        if checker.check_fg("Secondary alcohol", reactant) and checker.check_fg(
                            "Ketone", product
                        ):
                            print(f"Found secondary alcohol in reactant and ketone in product")
                            found_pattern = True
                            break

                elif reduction_reaction:
                    print("Found ketone reduction reaction")

                    # In forward direction: ketone → alcohol
                    # But in retrosynthesis, this means the final product was made by reducing a ketone
                    # which is equivalent to a late-stage oxidation strategy
                    for reactant in reactants:
                        if checker.check_fg("Ketone", reactant) and checker.check_fg(
                            "Secondary alcohol", product
                        ):
                            print(f"Found ketone in reactant and secondary alcohol in product")
                            found_pattern = True
                            break

                else:
                    print(
                        "Checking for functional group transformations without specific reaction type"
                    )
                    # Check for alcohol → ketone transformation
                    alcohol_to_ketone = False
                    for reactant in reactants:
                        if checker.check_fg("Secondary alcohol", reactant) and checker.check_fg(
                            "Ketone", product
                        ):
                            print(f"Found secondary alcohol in reactant and ketone in product")
                            alcohol_to_ketone = True
                            break

                    # Check for ketone → alcohol transformation
                    ketone_to_alcohol = False
                    for reactant in reactants:
                        if checker.check_fg("Ketone", reactant) and checker.check_fg(
                            "Secondary alcohol", product
                        ):
                            print(f"Found ketone in reactant and secondary alcohol in product")
                            ketone_to_alcohol = True
                            break

                    # In retrosynthesis, ketone → alcohol means the final product was made by reducing a ketone
                    if alcohol_to_ketone or ketone_to_alcohol:
                        found_pattern = True

        # Process children
        for child in node.get("children", []):
            if not found_pattern:  # Stop traversal if pattern already found
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {found_pattern}")
    return found_pattern
