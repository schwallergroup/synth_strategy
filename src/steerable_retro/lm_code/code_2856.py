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
    This function detects if the synthesis route involves reduction of a ketone to an alcohol.
    In retrosynthetic analysis, this corresponds to oxidation of an alcohol to a ketone.
    """
    ketone_reduction_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal ketone_reduction_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Primary check: Use the specific reaction checker
            # In retrosynthesis, check for both reduction and oxidation reactions
            if checker.check_reaction(
                "Reduction of aldehydes and ketones to alcohols", rsmi
            ) or checker.check_reaction(
                "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
            ):
                print(f"Ketone reduction/alcohol oxidation reaction detected: {rsmi}")
                ketone_reduction_detected = True
            else:
                # Secondary check: Look for alcohol to ketone transformation (retrosynthetic)
                # or ketone to alcohol transformation (forward)
                alcohol_in_product = (
                    checker.check_fg("Secondary alcohol", product)
                    or checker.check_fg("Primary alcohol", product)
                    or checker.check_fg("Tertiary alcohol", product)
                )
                ketone_in_product = checker.check_fg("Ketone", product)

                if alcohol_in_product:
                    print(f"Found alcohol in product: {product}")
                    for reactant in reactants:
                        if checker.check_fg("Ketone", reactant):
                            print(f"Found ketone in reactant: {reactant}")
                            ketone_reduction_detected = True
                            print(f"Ketone reduction detected through functional group analysis")
                            break

                # Also check the reverse direction (for completeness)
                if ketone_in_product and not ketone_reduction_detected:
                    print(f"Found ketone in product: {product}")
                    for reactant in reactants:
                        if (
                            checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Tertiary alcohol", reactant)
                        ):
                            print(f"Found alcohol in reactant: {reactant}")
                            ketone_reduction_detected = True
                            print(
                                f"Alcohol oxidation (reverse of ketone reduction) detected through functional group analysis"
                            )
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return ketone_reduction_detected
