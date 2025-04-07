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
    Detects if the route contains a late-stage C-N bond formation.
    """
    cn_bond_formation_found = False
    max_depth = 0

    # First pass: determine the maximum depth of the synthesis tree
    def calculate_max_depth(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)
        for child in node.get("children", []):
            calculate_max_depth(child, depth + 1)

    calculate_max_depth(route)
    print(f"Maximum synthesis depth: {max_depth}")

    # Second pass: check for C-N bond forming reactions
    def dfs_traverse(node, depth=0):
        nonlocal cn_bond_formation_found

        # Define what counts as "late stage" (first third of synthesis)
        is_late_stage = depth <= max_depth // 3

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Check for C-N bond forming reactions
            cn_bond_forming_reactions = [
                "N-alkylation of primary amines with alkyl halides",
                "N-alkylation of secondary amines with alkyl halides",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                "Reductive amination with aldehyde",
                "Reductive amination with ketone",
                "Reductive amination with alcohol",
                "Acylation of primary amines",
                "Acylation of secondary amines",
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                "Ugi reaction",
                "Goldberg coupling",
                "Ullmann-Goldberg Substitution amine",
                "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                "aza-Michael addition aromatic",
                "aza-Michael addition secondary",
                "aza-Michael addition primary",
            ]

            # Check for specific reaction types
            specific_reaction_found = False
            for reaction_type in cn_bond_forming_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Found C-N bond forming reaction: {reaction_type}")
                    specific_reaction_found = True

                    if is_late_stage:
                        print(f"Late-stage C-N bond formation at depth {depth}/{max_depth}")
                        cn_bond_formation_found = True
                        break

            # Always check for general C-N bond formation if we're in the late stage
            # and no specific reaction was identified
            if not cn_bond_formation_found and is_late_stage:
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for nitrogen-containing functional groups in reactants
                    has_nitrogen_fg = False
                    for reactant in reactants:
                        if reactant and (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Tertiary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                            or checker.check_fg("Azide", reactant)
                            or checker.check_fg("Amidinium", reactant)
                            or checker.check_fg("Primary amide", reactant)
                            or checker.check_fg("Secondary amide", reactant)
                            or checker.check_fg("Tertiary amide", reactant)
                        ):
                            has_nitrogen_fg = True
                            print(
                                f"Found nitrogen-containing functional group in reactant: {reactant}"
                            )
                            break

                    if has_nitrogen_fg:
                        # Check for new C-N bond formation by comparing reactants and product
                        prod_mol = Chem.MolFromSmiles(product)

                        # Define multiple SMARTS patterns for different types of C-N bonds
                        cn_bond_patterns = [
                            Chem.MolFromSmarts("[C]-[N]"),  # Single bond
                            Chem.MolFromSmarts("[C]=[N]"),  # Double bond
                            Chem.MolFromSmarts("[C]#[N]"),  # Triple bond
                            Chem.MolFromSmarts("[C]:[N]"),  # Aromatic bond
                        ]

                        # Count C-N bonds in product
                        prod_cn_count = 0
                        if prod_mol:
                            for pattern in cn_bond_patterns:
                                if prod_mol.HasSubstructMatch(pattern):
                                    prod_cn_count += len(prod_mol.GetSubstructMatches(pattern))

                        # Count C-N bonds in reactants
                        reactant_cn_count = 0
                        for reactant in reactants:
                            if reactant:
                                react_mol = Chem.MolFromSmiles(reactant)
                                if react_mol:
                                    for pattern in cn_bond_patterns:
                                        if react_mol.HasSubstructMatch(pattern):
                                            reactant_cn_count += len(
                                                react_mol.GetSubstructMatches(pattern)
                                            )

                        print(
                            f"C-N bonds in product: {prod_cn_count}, C-N bonds in reactants: {reactant_cn_count}"
                        )

                        # If product has more C-N bonds than reactants combined, a new bond was formed
                        if prod_cn_count > reactant_cn_count:
                            print(f"Found general C-N bond formation at depth {depth}")
                            cn_bond_formation_found = True
                except Exception as e:
                    print(f"Error in general C-N bond detection: {e}")

        for child in node.get("children", []):
            if (
                not cn_bond_formation_found
            ):  # Stop traversal if we already found what we're looking for
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if not cn_bond_formation_found:
        print("No late-stage C-N bond formation found in the route")

    return cn_bond_formation_found
