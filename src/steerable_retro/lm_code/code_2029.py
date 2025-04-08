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
    This function detects bidirectional alcohol-ketone interconversions.
    It looks for both alcohol→ketone and ketone→alcohol transformations in the same route.
    """
    alcohol_to_ketone = False
    ketone_to_alcohol = False

    print("Starting bidirectional alcohol-ketone interconversion check")

    def dfs_traverse(node):
        nonlocal alcohol_to_ketone, ketone_to_alcohol

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                print(f"Checking reaction: {rsmi}")

                try:
                    # Split reactants for individual checking
                    reactants = reactants_str.split(".")

                    # Check for alcohol to ketone/aldehyde/carboxylic acid/ester (oxidation)
                    if (
                        checker.check_reaction(
                            "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
                            rsmi,
                        )
                        or checker.check_reaction("Oxidation of alcohol to carboxylic acid", rsmi)
                        or checker.check_reaction(
                            "Oxidation of alcohol and aldehyde to ester", rsmi
                        )
                        or checker.check_reaction(
                            "Oxidative esterification of primary alcohols", rsmi
                        )
                        or checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    ):

                        # Verify that an alcohol is present in reactants
                        has_alcohol = any(
                            checker.check_fg("Primary alcohol", r)
                            or checker.check_fg("Secondary alcohol", r)
                            or checker.check_fg("Tertiary alcohol", r)
                            or checker.check_fg("Aromatic alcohol", r)
                            or checker.check_fg("Enol", r)
                            for r in reactants
                        )

                        # Verify that a carbonyl compound is present in product
                        has_carbonyl = (
                            checker.check_fg("Ketone", product_str)
                            or checker.check_fg("Aldehyde", product_str)
                            or checker.check_fg("Carboxylic acid", product_str)
                            or checker.check_fg("Ester", product_str)
                        )

                        if has_alcohol and has_carbonyl:
                            print("Alcohol to carbonyl transformation detected")
                            alcohol_to_ketone = True

                    # Check for ketone/aldehyde/ester to alcohol (reduction)
                    if (
                        checker.check_reaction(
                            "Reduction of aldehydes and ketones to alcohols", rsmi
                        )
                        or checker.check_reaction("Reduction of ketone to secondary alcohol", rsmi)
                        or checker.check_reaction("Grignard from aldehyde to alcohol", rsmi)
                        or checker.check_reaction("Grignard from ketone to alcohol", rsmi)
                        or checker.check_reaction(
                            "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                        )
                        or checker.check_reaction(
                            "Ester saponification (methyl deprotection)", rsmi
                        )
                        or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                    ):

                        # Verify that a carbonyl compound is present in reactants
                        has_carbonyl = any(
                            checker.check_fg("Ketone", r)
                            or checker.check_fg("Aldehyde", r)
                            or checker.check_fg("Formaldehyde", r)
                            or checker.check_fg("Ester", r)
                            or checker.check_fg("Carboxylic acid", r)
                            for r in reactants
                        )

                        # Verify that an alcohol is present in product
                        has_alcohol = (
                            checker.check_fg("Primary alcohol", product_str)
                            or checker.check_fg("Secondary alcohol", product_str)
                            or checker.check_fg("Tertiary alcohol", product_str)
                            or checker.check_fg("Aromatic alcohol", product_str)
                        )

                        if has_carbonyl and has_alcohol:
                            print("Carbonyl to alcohol transformation detected")
                            ketone_to_alcohol = True

                    # Direct check for alcohol to carbonyl transformation
                    if not alcohol_to_ketone:
                        has_alcohol_reactant = any(
                            checker.check_fg("Primary alcohol", r)
                            or checker.check_fg("Secondary alcohol", r)
                            or checker.check_fg("Tertiary alcohol", r)
                            or checker.check_fg("Aromatic alcohol", r)
                            or checker.check_fg("Enol", r)
                            for r in reactants
                        )

                        has_carbonyl_product = (
                            checker.check_fg("Ketone", product_str)
                            or checker.check_fg("Aldehyde", product_str)
                            or checker.check_fg("Carboxylic acid", product_str)
                            or checker.check_fg("Ester", product_str)
                        )

                        if has_alcohol_reactant and has_carbonyl_product:
                            print("Direct alcohol to carbonyl transformation detected")
                            alcohol_to_ketone = True

                    # Direct check for carbonyl to alcohol transformation
                    if not ketone_to_alcohol:
                        has_carbonyl_reactant = any(
                            checker.check_fg("Ketone", r)
                            or checker.check_fg("Aldehyde", r)
                            or checker.check_fg("Formaldehyde", r)
                            or checker.check_fg("Ester", r)
                            or checker.check_fg("Carboxylic acid", r)
                            for r in reactants
                        )

                        has_alcohol_product = (
                            checker.check_fg("Primary alcohol", product_str)
                            or checker.check_fg("Secondary alcohol", product_str)
                            or checker.check_fg("Tertiary alcohol", product_str)
                            or checker.check_fg("Aromatic alcohol", product_str)
                        )

                        if has_carbonyl_reactant and has_alcohol_product:
                            print("Direct carbonyl to alcohol transformation detected")
                            ketone_to_alcohol = True

                except Exception as e:
                    print(f"Error processing reaction: {e}")

                print(
                    f"Current status - alcohol_to_ketone: {alcohol_to_ketone}, ketone_to_alcohol: {ketone_to_alcohol}"
                )

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(
        f"Final result - alcohol_to_ketone: {alcohol_to_ketone}, ketone_to_alcohol: {ketone_to_alcohol}"
    )
    # Return True if both transformations are found
    return alcohol_to_ketone and ketone_to_alcohol
