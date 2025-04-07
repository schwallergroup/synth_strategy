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
    Detects if the synthetic route concludes with an esterification of a carboxylic acid
    or a saponification of an ester (which is the reverse of esterification)
    """
    final_step_is_esterification_or_saponification = False

    # First, find the depth of the first reaction node
    first_reaction_depth = None

    def find_first_reaction(node, depth=0):
        nonlocal first_reaction_depth
        if node["type"] == "reaction" and (
            first_reaction_depth is None or depth < first_reaction_depth
        ):
            first_reaction_depth = depth

        for child in node.get("children", []):
            find_first_reaction(child, depth + 1)

    find_first_reaction(route)
    print(f"First reaction depth: {first_reaction_depth}")

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_esterification_or_saponification

        # Check if this is a reaction node at the first level (final step in synthesis)
        if node["type"] == "reaction" and depth == first_reaction_depth:
            print(f"Examining final reaction step at depth {depth}")

            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Reaction SMILES: {rsmi}")

                # Check if this is an esterification or saponification reaction
                if checker.check_reaction("Esterification of Carboxylic Acids", rsmi):
                    print("Detected 'Esterification of Carboxylic Acids' reaction")
                    final_step_is_esterification_or_saponification = True
                elif checker.check_reaction("Transesterification", rsmi):
                    print("Detected 'Transesterification' reaction")
                    final_step_is_esterification_or_saponification = True
                elif checker.check_reaction("Ester saponification (methyl deprotection)", rsmi):
                    print("Detected 'Ester saponification (methyl deprotection)' reaction")
                    final_step_is_esterification_or_saponification = True
                elif checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi):
                    print("Detected 'Ester saponification (alkyl deprotection)' reaction")
                    final_step_is_esterification_or_saponification = True
                elif checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                ):
                    print(
                        "Detected 'Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters' reaction"
                    )
                    final_step_is_esterification_or_saponification = True
                else:
                    # Fallback method: check for functional group transformations in both directions
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Reactants: {reactants}")
                    print(f"Product: {product}")

                    # Check for esterification (carboxylic acid → ester)
                    has_carboxylic_acid_in_reactants = False
                    for reactant in reactants:
                        if checker.check_fg("Carboxylic acid", reactant):
                            print(f"Found carboxylic acid in reactant: {reactant}")
                            has_carboxylic_acid_in_reactants = True

                    has_ester_in_product = checker.check_fg("Ester", product)
                    if has_ester_in_product:
                        print(f"Found ester in product: {product}")

                    carboxylic_acid_in_product = checker.check_fg("Carboxylic acid", product)

                    # Check for saponification (ester → carboxylic acid)
                    has_ester_in_reactants = False
                    for reactant in reactants:
                        if checker.check_fg("Ester", reactant):
                            print(f"Found ester in reactant: {reactant}")
                            has_ester_in_reactants = True

                    has_carboxylic_acid_in_product = checker.check_fg("Carboxylic acid", product)
                    if has_carboxylic_acid_in_product:
                        print(f"Found carboxylic acid in product: {product}")

                    ester_in_product = checker.check_fg("Ester", product)

                    # Determine if either esterification or saponification occurred
                    if (
                        has_carboxylic_acid_in_reactants
                        and has_ester_in_product
                        and not carboxylic_acid_in_product
                    ) or (
                        has_ester_in_reactants
                        and has_carboxylic_acid_in_product
                        and not ester_in_product
                    ):
                        print(
                            "Detected esterification or saponification based on functional group analysis"
                        )
                        final_step_is_esterification_or_saponification = True
            except (KeyError, TypeError) as e:
                print(f"Error processing reaction metadata: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Final result: {final_step_is_esterification_or_saponification}")
    return final_step_is_esterification_or_saponification
