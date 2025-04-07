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
    This function detects if the synthetic route involves both oxidation and reduction steps.
    """
    oxidation_found = False
    reduction_found = False

    # Lists of oxidation and reduction reaction types
    oxidation_reactions = [
        "Oxidation of aldehydes to carboxylic acids",
        "Oxidation of ketone to carboxylic acid",
        "Oxidation of alcohol to carboxylic acid",
        "Oxidation of primary alcohols",
        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
        "Oxidation of alkene to aldehyde",
        "Oxidation of alkene to carboxylic acid",
        "Oxidation of alcohol and aldehyde to ester",
        "Oxidative esterification of primary alcohols",
        "Quinone formation",
        "Aromatic hydroxylation",
        "Sulfanyl to sulfinyl_peroxide",
        "Sulfanyl to sulfinyl_H2O2",
    ]

    reduction_reactions = [
        "Reduction of aldehydes and ketones to alcohols",
        "Reduction of ester to primary alcohol",
        "Reduction of ketone to secondary alcohol",
        "Reduction of carboxylic acid to primary alcohol",
        "Reduction of nitro groups to amines",
        "Reduction of nitrile to amine",
        "Reduction of primary amides to amines",
        "Reduction of secondary amides to amines",
        "Reduction of tertiary amides to amines",
        "Nef reaction (nitro to ketone)",
        "Azide to amine reduction (Staudinger)",
    ]

    # Functional group pairs for oxidation/reduction
    oxidation_fg_pairs = [
        ("Primary alcohol", "Aldehyde"),
        ("Primary alcohol", "Carboxylic acid"),
        ("Secondary alcohol", "Ketone"),
        ("Aldehyde", "Carboxylic acid"),
        ("Primary amine", "Nitrile"),
        ("Primary amine", "Nitro group"),
        ("Aliphatic thiol", "Sulfoxide"),
        ("Aliphatic thiol", "Sulfone"),
    ]

    reduction_fg_pairs = [
        ("Aldehyde", "Primary alcohol"),
        ("Ketone", "Secondary alcohol"),
        ("Carboxylic acid", "Primary alcohol"),
        ("Ester", "Primary alcohol"),
        ("Nitrile", "Primary amine"),
        ("Nitro group", "Primary amine"),
        ("Azide", "Primary amine"),
        ("Sulfoxide", "Aliphatic thiol"),
        ("Sulfone", "Aliphatic thiol"),
    ]

    def dfs_traverse(node):
        nonlocal oxidation_found, reduction_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Check for known oxidation reactions
                for rxn_type in oxidation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Oxidation detected: {rxn_type} in reaction: {rsmi}")
                        oxidation_found = True
                        break

                # Check for known reduction reactions
                for rxn_type in reduction_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Reduction detected: {rxn_type} in reaction: {rsmi}")
                        reduction_found = True
                        break

                # If no specific reaction type was found, check for functional group transformations
                if not (oxidation_found and reduction_found):
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for oxidation based on functional group changes
                    if not oxidation_found:
                        for r_fg, p_fg in oxidation_fg_pairs:
                            if any(
                                checker.check_fg(r_fg, r) for r in reactants
                            ) and checker.check_fg(p_fg, product):
                                print(f"Oxidation detected: {r_fg} to {p_fg} in reaction: {rsmi}")
                                oxidation_found = True
                                break

                    # Check for reduction based on functional group changes
                    if not reduction_found:
                        for r_fg, p_fg in reduction_fg_pairs:
                            if any(
                                checker.check_fg(r_fg, r) for r in reactants
                            ) and checker.check_fg(p_fg, product):
                                print(f"Reduction detected: {r_fg} to {p_fg} in reaction: {rsmi}")
                                reduction_found = True
                                break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    result = oxidation_found and reduction_found
    print(f"Bidirectional redox strategy detected: {result}")
    return result
