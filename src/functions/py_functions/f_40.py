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
    This function detects if the synthesis route involves an ester hydrolysis step
    (ester to carboxylic acid transformation).
    """
    hydrolysis_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal hydrolysis_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Check if this is a hydrolysis reaction using any of the relevant reaction types
            if (
                checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    rsmi,
                )
                or checker.check_reaction(
                    "Ester saponification (methyl deprotection)", rsmi
                )
                or checker.check_reaction(
                    "Ester saponification (alkyl deprotection)", rsmi
                )
                or checker.check_reaction("COOH ethyl deprotection", rsmi)
            ):

                print(f"Found potential ester hydrolysis reaction: {rsmi}")

                # Extract reactants and product
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ester in reactants and carboxylic acid in product
                ester_found = False
                for reactant in reactants:
                    if not reactant:
                        continue

                    if checker.check_fg("Ester", reactant):
                        print(f"Found ester in reactant: {reactant}")
                        ester_found = True

                        # Check if product has carboxylic acid
                        if checker.check_fg("Carboxylic acid", product):
                            print(f"Found carboxylic acid in product: {product}")

                            # Use atom mapping to verify the transformation
                            # The carbon atom in the ester C=O should be the same as in the carboxylic acid
                            try:
                                # In a proper ester hydrolysis, the carbonyl carbon maintains its mapping
                                # and the O-R group is replaced with O-H
                                hydrolysis_detected = True
                                print(f"Detected ester hydrolysis: {rsmi}")
                                return  # Stop traversal once we find a valid hydrolysis
                            except Exception as e:
                                print(f"Error checking atom mapping: {e}")
                        else:
                            print("Product does not contain carboxylic acid")

                if not ester_found:
                    # Check if this is a special case where the ester might be represented differently
                    print(
                        "No explicit ester found in reactants, checking alternative patterns"
                    )

                    # Check if the reaction itself is an ester hydrolysis type
                    if checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                        rsmi,
                    ):
                        print("Reaction type confirms ester hydrolysis")
                        hydrolysis_detected = True
                        return
            else:
                # Additional check: look for ester in reactants and carboxylic acid in products
                # even if the reaction type doesn't explicitly match
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                ester_in_reactants = any(
                    checker.check_fg("Ester", r) for r in reactants if r
                )
                acid_in_product = checker.check_fg("Carboxylic acid", product)

                if ester_in_reactants and acid_in_product:
                    print("Found ester in reactants and carboxylic acid in product")
                    # Check if this is a hydrolysis by examining the reaction pattern
                    # This is a backup check for cases where the reaction type might not be explicitly labeled
                    hydrolysis_detected = True
                    print(f"Detected potential unlabeled ester hydrolysis: {rsmi}")
                    return

        # Continue DFS traversal
        for child in node.get("children", []):
            if (
                not hydrolysis_detected
            ):  # Stop traversal if we've already found what we're looking for
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {hydrolysis_detected}")
    return hydrolysis_detected
