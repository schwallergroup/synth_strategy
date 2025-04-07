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


def main(route, min_count=3):
    """
    This function detects if the route contains multiple C-N bond formations.
    """
    cn_bond_formation_count = 0

    # Comprehensive list of C-N bond formation reaction types
    cn_formation_reactions = [
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Acyl chloride with ammonia to amide",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with secondary amine to amide",
        "Carboxylic acid with primary amine to amide",
        "Ester with ammonia to amide",
        "Ester with primary amine to amide",
        "Ester with secondary amine to amide",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "Reductive amination with alcohol",
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Alkylation of amines",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Ugi reaction",
        "Mignonac reaction",
        "Schmidt reaction acid_amine",
        "Schmidt ketone amide",
        "Schmidt aldehyde amide",
        "aza-Michael addition aromatic",
        "aza-Michael addition secondary",
        "aza-Michael addition primary",
        "Acylation of secondary amines",
        "Acylation of primary amines",
        "Displacement of ethoxy group by primary amine",
        "Displacement of ethoxy group by secondary amine",
        "Urea synthesis via isocyanate and primary amine",
        "Urea synthesis via isocyanate and secondary amine",
        "Urea synthesis via isocyanate and diazo",
        "Urea synthesis via isocyanate and sulfonamide",
        "{Buchwald-Hartwig}",
        "Sulfonamide synthesis (Schotten-Baumann) primary amine",
        "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
        "Schotten-Baumann to ester",
        "{Schotten-Baumann_amide}",
        "{sulfon_amide}",
        "{urea}",
        "{N-arylation_heterocycles}",
    ]

    def dfs_traverse(node):
        nonlocal cn_bond_formation_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if this reaction is a C-N bond formation using the checker function
            is_cn_formation = False
            for rxn_type in cn_formation_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    is_cn_formation = True
                    print(f"C-N bond formation detected: {rxn_type} in reaction: {rsmi}")
                    break

            # If not detected by reaction type, check for specific functional group changes
            if not is_cn_formation:
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for amide formation
                product_has_amide = (
                    checker.check_fg("Primary amide", product_smiles)
                    or checker.check_fg("Secondary amide", product_smiles)
                    or checker.check_fg("Tertiary amide", product_smiles)
                )

                reactants_have_amide = any(
                    checker.check_fg("Primary amide", r)
                    or checker.check_fg("Secondary amide", r)
                    or checker.check_fg("Tertiary amide", r)
                    for r in reactants_smiles
                )

                # Check for amine formation
                product_has_amine = (
                    checker.check_fg("Primary amine", product_smiles)
                    or checker.check_fg("Secondary amine", product_smiles)
                    or checker.check_fg("Tertiary amine", product_smiles)
                )

                reactants_have_amine = any(
                    checker.check_fg("Primary amine", r)
                    or checker.check_fg("Secondary amine", r)
                    or checker.check_fg("Tertiary amine", r)
                    for r in reactants_smiles
                )

                # Check for other nitrogen-containing functional groups
                product_has_urea = checker.check_fg("Urea", product_smiles)
                reactants_have_urea = any(checker.check_fg("Urea", r) for r in reactants_smiles)

                product_has_sulfonamide = checker.check_fg("Sulfonamide", product_smiles)
                reactants_have_sulfonamide = any(
                    checker.check_fg("Sulfonamide", r) for r in reactants_smiles
                )

                if (
                    (product_has_amide and not reactants_have_amide)
                    or (product_has_amine and not reactants_have_amine)
                    or (product_has_urea and not reactants_have_urea)
                    or (product_has_sulfonamide and not reactants_have_sulfonamide)
                ):
                    is_cn_formation = True
                    print(
                        f"C-N bond formation detected through functional group analysis in reaction: {rsmi}"
                    )

            if is_cn_formation:
                cn_bond_formation_count += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"C-N bond formation count: {cn_bond_formation_count}")
    return cn_bond_formation_count >= min_count
