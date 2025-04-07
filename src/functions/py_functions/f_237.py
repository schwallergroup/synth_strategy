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
    Detects if the synthesis route contains an aryl halide to ester conversion.
    """
    conversion_found = False

    # Track molecules with aryl halides for multi-step conversion detection
    aryl_halide_molecules = {}

    def dfs_traverse(node, depth=0, path=None):
        nonlocal conversion_found
        if path is None:
            path = []

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Track molecules with aryl halides
            if checker.check_fg("Aromatic halide", mol_smiles):
                aryl_halide_molecules[mol_smiles] = depth
                print(f"Found molecule with aryl halide at depth {depth}: {mol_smiles}")

            # Check if this molecule has an ester and was derived from an aryl halide
            if checker.check_fg("Ester", mol_smiles) and path:
                print(f"Found molecule with ester at depth {depth}: {mol_smiles}")
                # Check if any parent in the path had an aryl halide
                for parent_smiles in path:
                    if parent_smiles in aryl_halide_molecules:
                        print(
                            f"Multi-step conversion detected: aryl halide at depth {aryl_halide_molecules[parent_smiles]} -> ester at depth {depth}"
                        )
                        conversion_found = True

        elif node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if any reactant contains an aromatic halide
            aryl_halide_reactants = [
                r for r in reactants if checker.check_fg("Aromatic halide", r)
            ]

            # Check if product contains an ester
            has_ester_product = checker.check_fg("Ester", product)

            if aryl_halide_reactants and has_ester_product:
                print(
                    f"Found potential aryl halide to ester conversion at depth {depth}"
                )

                # Expanded list of relevant reactions
                relevant_reactions = [
                    "Carbonylation with aryl formates",
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic acids OTf",
                    "Suzuki coupling with boronic esters",
                    "Suzuki coupling with boronic esters OTf",
                    "Suzuki coupling with sulfonic esters",
                    "Oxidative esterification of primary alcohols",
                    "Esterification of Carboxylic Acids",
                    "Heck reaction with vinyl ester",
                    "Heck reaction with vinyl ester and amine",
                    "Oxidative Heck reaction with vinyl ester",
                    "Sonogashira alkyne_aryl halide",
                    "Sonogashira acetylene_aryl halide",
                    "Sonogashira alkyne_aryl OTf",
                    "Sonogashira acetylene_aryl OTf",
                    "Schotten-Baumann to ester",
                    "Transesterification",
                    "O-alkylation of carboxylic acids with diazo compounds",
                ]

                # Check for any relevant reaction type
                for reaction_type in relevant_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Matched reaction type: {reaction_type}")
                        conversion_found = True
                        break

                # If no specific reaction type matched, try to verify the transformation
                if not conversion_found:
                    print(
                        "No specific reaction type matched, checking general transformation..."
                    )
                    # Check if the reaction involves carbonylation or esterification
                    # This is a fallback for reactions not in our list
                    for r in aryl_halide_reactants:
                        if not checker.check_fg("Ester", r) and has_ester_product:
                            print(
                                "General aryl halide to ester transformation detected"
                            )
                            conversion_found = True
                            break

        # Add current molecule to path for multi-step tracking
        if node["type"] == "mol":
            new_path = path + [node["smiles"]]
        else:
            new_path = path

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, new_path)

    dfs_traverse(route)
    return conversion_found
