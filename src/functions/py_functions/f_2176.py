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
    Detects if the synthesis includes esterification reactions
    (carboxylic acid → ester transformation).
    """
    esterification_detected = False

    def dfs_traverse(node):
        nonlocal esterification_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction: {rsmi}")

            # Check for various esterification reaction types
            esterification_reactions = [
                "Esterification of Carboxylic Acids",
                "Transesterification",
                "O-alkylation of carboxylic acids with diazo compounds",
                "Oxidative esterification of primary alcohols",
                "Acetic anhydride and alcohol to ester",
                "Pinner reaction to ester",
                "Schotten-Baumann to ester",
            ]

            # Check if any of the esterification reaction types match
            if any(
                checker.check_reaction(reaction_type, rsmi)
                for reaction_type in esterification_reactions
            ):
                print(f"Detected esterification reaction: {rsmi}")
                esterification_detected = True

            # If not identified by reaction type, check by functional group transformation
            try:
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
                product_mol = Chem.MolFromSmiles(product_part)

                if product_mol and all(r is not None for r in reactants):
                    # Check for carboxylic acid in reactants
                    carboxylic_reactants = []
                    for i, r in enumerate(reactants):
                        r_smiles = Chem.MolToSmiles(r)
                        if checker.check_fg("Carboxylic acid", r_smiles):
                            print(f"Found carboxylic acid in reactant {i}: {r_smiles}")
                            carboxylic_reactants.append(r_smiles)

                    # Check for ester in reactants (for transesterification)
                    ester_reactants = []
                    for i, r in enumerate(reactants):
                        r_smiles = Chem.MolToSmiles(r)
                        if checker.check_fg("Ester", r_smiles):
                            print(f"Found ester in reactant {i}: {r_smiles}")
                            ester_reactants.append(r_smiles)

                    # Check for alcohol in reactants
                    alcohol_reactants = []
                    for i, r in enumerate(reactants):
                        r_smiles = Chem.MolToSmiles(r)
                        if any(
                            checker.check_fg(fg, r_smiles)
                            for fg in [
                                "Primary alcohol",
                                "Secondary alcohol",
                                "Tertiary alcohol",
                                "Aromatic alcohol",
                            ]
                        ):
                            print(f"Found alcohol in reactant {i}: {r_smiles}")
                            alcohol_reactants.append(r_smiles)

                    # Check for ester in product
                    product_smiles = Chem.MolToSmiles(product_mol)
                    has_ester_product = checker.check_fg("Ester", product_smiles)
                    if has_ester_product:
                        print(f"Found ester in product: {product_smiles}")

                    # Check for carboxylic acid in product (for ester hydrolysis)
                    has_carboxylic_product = checker.check_fg(
                        "Carboxylic acid", product_smiles
                    )
                    if has_carboxylic_product:
                        print(f"Found carboxylic acid in product: {product_smiles}")

                    # Case 1: Carboxylic acid to ester (direct esterification)
                    if carboxylic_reactants and has_ester_product:
                        print(
                            f"Detected esterification (carboxylic acid to ester): {rsmi}"
                        )
                        esterification_detected = True

                    # Case 2: Alcohol + carboxylic acid to ester
                    if alcohol_reactants and carboxylic_reactants and has_ester_product:
                        print(
                            f"Detected alcohol + carboxylic acid esterification: {rsmi}"
                        )
                        esterification_detected = True

                    # Case 3: Ester to different ester (transesterification)
                    if ester_reactants and has_ester_product and alcohol_reactants:
                        print(f"Detected transesterification (ester to ester): {rsmi}")
                        esterification_detected = True

                    # Case 4: Ester hydrolysis (reverse of esterification)
                    if ester_reactants and has_carboxylic_product:
                        print(
                            f"Detected ester hydrolysis (ester to carboxylic acid): {rsmi}"
                        )
                        esterification_detected = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Esterification detected: {esterification_detected}")
    return esterification_detected
