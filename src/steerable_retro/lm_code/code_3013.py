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
    This function detects late-stage mesylation (in the final step of synthesis).
    Looks for transformation of alcohol to mesylate in the first reaction step.
    """
    mesylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal mesylation_detected

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        # Check if this is a reaction node at the late stage (depth 0 or 1)
        if node["type"] == "reaction" and (
            depth == 1 or (depth == 0 and "rsmi" in node.get("metadata", {}))
        ):
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for alcohol in reactants
                has_alcohol = False
                alcohol_reactant = None
                for reactant in reactants:
                    if (
                        checker.check_fg("Primary alcohol", reactant)
                        or checker.check_fg("Secondary alcohol", reactant)
                        or checker.check_fg("Tertiary alcohol", reactant)
                        or checker.check_fg("Aromatic alcohol", reactant)
                    ):
                        has_alcohol = True
                        alcohol_reactant = reactant
                        print(f"Found alcohol in reactant: {reactant}")
                        break

                # Check for mesylate in product
                has_mesylate = checker.check_fg("Mesylate", product)
                if has_mesylate:
                    print(f"Found mesylate in product: {product}")
                else:
                    print(f"No mesylate found in product: {product}")

                # Check if this is a mesylation reaction using multiple methods
                is_mesylation_reaction = checker.check_reaction(
                    "Formation of Sulfonic Esters", rsmi
                ) or checker.check_reaction(
                    "Formation of Sulfonic Esters on TMS protected alcohol", rsmi
                )

                if is_mesylation_reaction:
                    print("Confirmed mesylation reaction type")
                else:
                    print("Not a mesylation reaction according to reaction checker")

                    # Manual check for mesylation reagents if reaction check fails
                    has_mesylating_agent = False
                    for reactant in reactants:
                        # Check for methanesulfonyl chloride or similar reagents
                        if (
                            "S(=O)(=O)Cl" in reactant
                            or "ClS(=O)(=O)" in reactant
                            or "S(=O)(=O)(Cl)" in reactant
                            or "S(C)(=O)(=O)Cl" in reactant
                            or checker.check_fg("Sulfonyl halide", reactant)
                        ):
                            has_mesylating_agent = True
                            print(f"Found mesylating agent: {reactant}")
                            break

                    # If we have alcohol, mesylate, and a mesylating agent, it's likely a mesylation
                    if has_alcohol and has_mesylate and has_mesylating_agent:
                        is_mesylation_reaction = True
                        print("Manually confirmed as mesylation reaction")

                # Final determination
                if (
                    has_alcohol
                    and has_mesylate
                    and (is_mesylation_reaction or "S(C)(=O)(=O)" in product)
                ):
                    print("Detected late-stage mesylation")
                    mesylation_detected = True
                else:
                    if not has_alcohol:
                        print("No alcohol found in reactants")
                    if not has_mesylate:
                        print("No mesylate found in product")
                    if not is_mesylation_reaction:
                        print("Not a mesylation reaction")

        # Continue DFS traversal
        for child in node.get("children", []):
            if not mesylation_detected:  # Stop traversal if we already found what we're looking for
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return mesylation_detected
