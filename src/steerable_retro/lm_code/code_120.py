#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    Detects if the synthesis route contains both oxidation and reduction steps.
    """
    found_oxidation = False
    found_reduction = False

    # List of oxidation reaction types
    oxidation_reactions = [
        "Oxidation of aldehydes to carboxylic acids",
        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
        "Oxidation of ketone to carboxylic acid",
        "Oxidation of alcohol to carboxylic acid",
        "Oxidation of nitrile to carboxylic acid",
        "Oxidation of amide to carboxylic acid",
        "Alkene oxidation to aldehyde",
        "Oxidative esterification of primary alcohols",
        "Oxidation of alcohol and aldehyde to ester",
        "Aromatic hydroxylation",
        "Oxidation of boronic acids",
        "Oxidation of boronic esters",
    ]

    # List of reduction reaction types
    reduction_reactions = [
        "Reduction of aldehydes and ketones to alcohols",
        "Reduction of ester to primary alcohol",
        "Reduction of ketone to secondary alcohol",
        "Reduction of carboxylic acid to primary alcohol",
        "Reduction of nitro groups to amines",
        "Reduction of nitrile to amide",
        "Reduction of nitrile to amine",
        "Reduction of primary amides to amines",
        "Reduction of secondary amides to amines",
        "Reduction of tertiary amides to amines",
        "Azide to amine reduction (Staudinger)",
    ]

    def dfs_traverse(node):
        nonlocal found_oxidation, found_reduction

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]

            try:
                # Check for oxidation reactions
                for rxn_type in oxidation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        found_oxidation = True
                        print(f"Found oxidation step: {rxn_type}")
                        break

                # Check for reduction reactions
                for rxn_type in reduction_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        found_reduction = True
                        print(f"Found reduction step: {rxn_type}")
                        break

                # If we haven't found specific reaction types, check for functional group changes
                if not found_oxidation or not found_reduction:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for oxidation by functional group changes
                    if not found_oxidation:
                        # Alcohol to aldehyde/ketone/acid/ester (oxidation)
                        alcohol_in_reactants = any(
                            checker.check_fg("Primary alcohol", r)
                            or checker.check_fg("Secondary alcohol", r)
                            or checker.check_fg("Tertiary alcohol", r)
                            or checker.check_fg("Aromatic alcohol", r)
                            for r in reactants
                        )

                        oxidized_in_product = (
                            checker.check_fg("Aldehyde", product)
                            or checker.check_fg("Ketone", product)
                            or checker.check_fg("Carboxylic acid", product)
                            or checker.check_fg("Ester", product)
                        )

                        if alcohol_in_reactants and oxidized_in_product:
                            found_oxidation = True
                            print("Found oxidation step (alcohol to carbonyl)")

                        # Aldehyde to carboxylic acid (oxidation)
                        aldehyde_in_reactants = any(
                            checker.check_fg("Aldehyde", r) for r in reactants
                        )
                        acid_in_product = checker.check_fg("Carboxylic acid", product)

                        if aldehyde_in_reactants and acid_in_product:
                            found_oxidation = True
                            print("Found oxidation step (aldehyde to acid)")

                    # Check for reduction by functional group changes
                    if not found_reduction:
                        # Carbonyl to alcohol (reduction)
                        carbonyl_in_reactants = any(
                            checker.check_fg("Aldehyde", r)
                            or checker.check_fg("Ketone", r)
                            or checker.check_fg("Carboxylic acid", r)
                            or checker.check_fg("Ester", r)
                            or checker.check_fg("Amide", r)
                            for r in reactants
                        )

                        alcohol_in_product = (
                            checker.check_fg("Primary alcohol", product)
                            or checker.check_fg("Secondary alcohol", product)
                            or checker.check_fg("Tertiary alcohol", product)
                        )

                        if carbonyl_in_reactants and alcohol_in_product:
                            found_reduction = True
                            print("Found reduction step (carbonyl to alcohol)")

                        # Nitro to amine (reduction)
                        nitro_in_reactants = any(
                            checker.check_fg("Nitro group", r) for r in reactants
                        )
                        amine_in_product = (
                            checker.check_fg("Primary amine", product)
                            or checker.check_fg("Secondary amine", product)
                            or checker.check_fg("Tertiary amine", product)
                            or checker.check_fg("Aniline", product)
                        )

                        if nitro_in_reactants and amine_in_product:
                            found_reduction = True
                            print("Found reduction step (nitro to amine)")

                        # Azide to amine (reduction)
                        azide_in_reactants = any(checker.check_fg("Azide", r) for r in reactants)

                        if azide_in_reactants and amine_in_product:
                            found_reduction = True
                            print("Found reduction step (azide to amine)")

            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we found both oxidation and reduction
    strategy_present = found_oxidation and found_reduction

    print(f"Oxidation-reduction sequence strategy detected: {strategy_present}")
    return strategy_present
