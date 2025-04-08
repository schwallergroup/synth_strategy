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
    This function detects if there is a C-N bond formation in the late stage of synthesis (depth <= 1).
    This includes both new C-N bond formation and C-N bond modifications (like nitro reduction).
    """
    has_late_stage_cn_bond = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_cn_bond

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            # Check if this is a late stage reaction (depth 0 or 1)
            depth_info = node.get("metadata", {}).get("ID", "")
            is_late_stage = False

            # Try to extract depth from ID field
            if any(f"Depth: {i}" in depth_info for i in range(2)):
                is_late_stage = True
            # Fallback: use the traversal depth
            elif depth <= 1:
                is_late_stage = True

            if is_late_stage:
                rsmi = node["metadata"]["rsmi"]
                print(f"Analyzing late-stage reaction: {rsmi}")

                # Check for known C-N bond formation or modification reactions
                cn_bond_reactions = [
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                    "Buchwald-Hartwig",
                    "Ullmann-Goldberg Substitution amine",
                    "Ullmann-Goldberg coupling",
                    "Goldberg coupling",
                    "Goldberg coupling aryl amine-aryl chloride",
                    "Goldberg coupling aryl amide-aryl chloride",
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "Reductive amination with alcohol",
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Ugi reaction",
                    "Schotten-Baumann_amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Aminolysis of esters",
                    "aza-Michael addition primary",
                    "aza-Michael addition secondary",
                    "aza-Michael addition aromatic",
                    "Petasis reaction with amines and boronic acids",
                    "Petasis reaction with amines and boronic esters",
                    "Petasis reaction with amines aldehydes and boronic acids",
                    "Reduction of nitro groups to amines",  # Added nitro reduction
                ]

                for reaction_type in cn_bond_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Detected late-stage C-N bond formation/modification: {reaction_type}"
                        )
                        has_late_stage_cn_bond = True
                        return

                # If no known reaction type is found, check for C-N bond formation or modification manually
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for specific functional group transformations
                # Check for nitro to amine transformation
                has_nitro_in_reactants = any(checker.check_fg("Nitro group", r) for r in reactants)
                has_amine_in_product = checker.check_fg(
                    "Primary amine", product
                ) or checker.check_fg("Aniline", product)

                print(
                    f"Manual check - Nitro in reactants: {has_nitro_in_reactants}, Amine in product: {has_amine_in_product}"
                )

                if has_nitro_in_reactants and has_amine_in_product:
                    print("Detected late-stage nitro to amine transformation")
                    has_late_stage_cn_bond = True
                    return

                # Check for reactants with potential C-N bond formation
                carbon_electrophiles = [
                    "Aromatic halide",
                    "Primary halide",
                    "Secondary halide",
                    "Tertiary halide",
                    "Acyl halide",
                    "Aldehyde",
                    "Ketone",
                ]

                nitrogen_nucleophiles = [
                    "Primary amine",
                    "Secondary amine",
                    "Tertiary amine",
                    "Aniline",
                    "Amide",
                    "Primary amide",
                    "Secondary amide",
                    "Tertiary amide",
                    "Azide",
                    "Hydrazine",
                ]

                has_carbon_electrophile = False
                has_nitrogen_nucleophile = False

                for r in reactants:
                    for fg in carbon_electrophiles:
                        if checker.check_fg(fg, r):
                            print(f"Found carbon electrophile: {fg} in {r}")
                            has_carbon_electrophile = True

                    for fg in nitrogen_nucleophiles:
                        if checker.check_fg(fg, r):
                            print(f"Found nitrogen nucleophile: {fg} in {r}")
                            has_nitrogen_nucleophile = True

                if has_carbon_electrophile and has_nitrogen_nucleophile:
                    # Verify C-N bond formation by checking product
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
                    ]

                    if product_mol:
                        # Check for various C-N bond patterns
                        cn_patterns = [
                            Chem.MolFromSmarts("C-N"),  # Aliphatic C-N
                            Chem.MolFromSmarts("c-N"),  # Aromatic C-aliphatic N
                            Chem.MolFromSmarts("C-n"),  # Aliphatic C-aromatic N
                            Chem.MolFromSmarts("c-n"),  # Aromatic C-aromatic N
                            Chem.MolFromSmarts("C=N"),  # C=N bond
                            Chem.MolFromSmarts("C#N"),  # C≡N bond (excluding nitriles)
                        ]

                        for pattern in cn_patterns:
                            if product_mol.HasSubstructMatch(pattern):
                                # Check if this is a new bond not present in reactants
                                new_bond = True
                                for mol in reactant_mols:
                                    if mol and mol.HasSubstructMatch(pattern):
                                        matches_in_product = len(
                                            product_mol.GetSubstructMatches(pattern)
                                        )
                                        matches_in_reactant = len(mol.GetSubstructMatches(pattern))
                                        if matches_in_product <= matches_in_reactant:
                                            new_bond = False
                                            break

                                if new_bond:
                                    print(
                                        f"Detected late-stage C-N bond formation via manual check"
                                    )
                                    has_late_stage_cn_bond = True
                                    return

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return has_late_stage_cn_bond
