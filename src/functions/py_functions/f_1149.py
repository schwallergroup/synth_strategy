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
    This function detects if esterification occurs at the final stages (depths 0-2) of the synthesis.
    In retrosynthetic analysis, esterification appears as ester hydrolysis in the forward direction.
    """
    late_stage_esterification_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_esterification_detected

        # Check at depths 0-2 (final reactions)
        if node["type"] == "reaction" and depth <= 2:
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is an esterification reaction (in retrosynthetic direction)
                is_esterification = False

                # Direct esterification reactions (forward direction)
                esterification_reactions = [
                    "Esterification of Carboxylic Acids",
                    "Transesterification",
                    "O-alkylation of carboxylic acids with diazo compounds",
                    "Oxidative esterification of primary alcohols",
                    "Acetic anhydride and alcohol to ester",
                    "Schotten-Baumann to ester",
                ]

                # Ester hydrolysis reactions (forward direction, which are esterifications in retrosynthesis)
                hydrolysis_reactions = [
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    "Ester saponification (methyl deprotection)",
                    "Ester saponification (alkyl deprotection)",
                    "COOH ethyl deprotection",
                ]

                # Check for direct esterification (forward direction)
                for reaction_type in esterification_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Detected forward esterification reaction type: {reaction_type}"
                        )
                        is_esterification = True
                        break

                # Check for ester hydrolysis (which is esterification in retrosynthesis)
                if not is_esterification:
                    for reaction_type in hydrolysis_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(
                                f"Detected ester hydrolysis reaction type: {reaction_type}"
                            )
                            is_esterification = True
                            break

                # If not directly identified, check for functional group changes
                if not is_esterification:
                    # Check for ester in reactants (retrosynthetic esterification)
                    ester_in_reactants = False
                    for reactant in reactants:
                        if checker.check_fg("Ester", reactant):
                            ester_in_reactants = True
                            print(f"Found ester in reactant: {reactant}")
                            break

                    # Check for carboxylic acid in product (retrosynthetic esterification)
                    acid_in_product = checker.check_fg("Carboxylic acid", product)
                    if acid_in_product:
                        print(f"Found carboxylic acid in product: {product}")

                    # Check for alcohol in product (might be present in hydrolysis)
                    alcohol_in_product = (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                        or checker.check_fg("Aromatic alcohol", product)
                    )

                    # Forward esterification check
                    acid_in_reactants = False
                    alcohol_in_reactants = False
                    for reactant in reactants:
                        if checker.check_fg("Carboxylic acid", reactant):
                            acid_in_reactants = True
                            print(f"Found carboxylic acid in reactant: {reactant}")
                        if (
                            checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Tertiary alcohol", reactant)
                            or checker.check_fg("Aromatic alcohol", reactant)
                        ):
                            alcohol_in_reactants = True
                            print(f"Found alcohol in reactant: {reactant}")

                    ester_in_product = checker.check_fg("Ester", product)
                    if ester_in_product:
                        print(f"Found ester in product: {product}")

                    # Determine if this is an esterification (in either direction)
                    if (ester_in_reactants and acid_in_product) or (
                        acid_in_reactants and alcohol_in_reactants and ester_in_product
                    ):
                        print(f"Detected esterification by functional group analysis")
                        is_esterification = True

                if is_esterification:
                    print(f"Late-stage esterification detected at depth {depth}")
                    late_stage_esterification_detected = True

        # Continue traversing the synthesis route
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Late-stage esterification detected: {late_stage_esterification_detected}")
    return late_stage_esterification_detected
