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
    Detects a strategy involving introduction and manipulation of a formyl group.
    """
    found_formyl_introduction = False
    found_formyl_manipulation = False

    # Store reactions for analysis
    formyl_introduction_reactions = []
    formyl_manipulation_reactions = []

    def dfs_traverse(node, depth=0):
        nonlocal found_formyl_introduction, found_formyl_manipulation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for formyl group presence
                product_has_formyl = checker.check_fg("Formaldehyde", product) or checker.check_fg(
                    "Aldehyde", product
                )
                reactants_have_formyl = any(
                    checker.check_fg("Formaldehyde", r) or checker.check_fg("Aldehyde", r)
                    for r in reactants
                )

                # Check for formyl introduction (typically early in synthesis)
                if depth >= 2:  # More flexible depth threshold
                    # Check for formylation reactions
                    is_formylation = (
                        checker.check_reaction("Carbonylation with aryl formates", rsmi)
                        or checker.check_reaction("Bouveault aldehyde synthesis", rsmi)
                        or checker.check_reaction("Oxidation of alkene to aldehyde", rsmi)
                        or checker.check_reaction("Oxidation of primary alcohol to aldehyde", rsmi)
                        or checker.check_reaction(
                            "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
                            rsmi,
                        )
                        or checker.check_reaction("Hydration of alkyne to aldehyde", rsmi)
                    )

                    # Check for common formylating agents in reactants
                    has_formylating_agent = any(
                        checker.check_fg("Formaldehyde", r)
                        or "CN(C)C=O" in r
                        or "OC=O" in r  # DMF
                        or "C(=O)Cl" in r  # Formic acid  # Formyl chloride
                        for r in reactants
                    )

                    if (
                        (product_has_formyl and not reactants_have_formyl)
                        or is_formylation
                        or (has_formylating_agent and product_has_formyl)
                    ):
                        found_formyl_introduction = True
                        formyl_introduction_reactions.append((rsmi, depth))
                        print(f"Found formyl group introduction at depth {depth}")
                        print(f"Reaction SMILES: {rsmi}")

                # Check for formyl group manipulation (can occur at any stage)
                # Check if formyl group is being transformed

                # Check for specific formyl manipulation reactions
                is_formyl_manipulation = (
                    checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rsmi)
                    or checker.check_reaction(
                        "Reduction of aldehydes and ketones to alcohols", rsmi
                    )
                    or checker.check_reaction("Aldol condensation", rsmi)
                    or checker.check_reaction("Wittig reaction with triphenylphosphorane", rsmi)
                    or checker.check_reaction("Reductive amination with aldehyde", rsmi)
                    or checker.check_reaction(
                        "Addition of primary amines to aldehydes/thiocarbonyls", rsmi
                    )
                    or checker.check_reaction(
                        "Addition of secondary amines to aldehydes/thiocarbonyls", rsmi
                    )
                    or checker.check_reaction("Benzothiazole formation from aldehyde", rsmi)
                    or checker.check_reaction("Benzoxazole formation from aldehyde", rsmi)
                    or checker.check_reaction("Benzimidazole formation from aldehyde", rsmi)
                    or checker.check_reaction("Homologation of aldehydes with formaldehyde", rsmi)
                    or checker.check_reaction(
                        "Derived alcohol via homologation of aldehydes with formaldehyde", rsmi
                    )
                    or checker.check_reaction(
                        "Derived benzimidazole via homologation of aldehydes with formaldehyde",
                        rsmi,
                    )
                    or checker.check_reaction("Acetal hydrolysis to aldehyde", rsmi)
                    or checker.check_reaction("Aldehyde or ketone acetalization", rsmi)
                )

                # Check if formyl group is being used in a reaction (present in reactants)
                formyl_used_in_reaction = reactants_have_formyl and (
                    # Check if the reaction is using the formyl group
                    not product_has_formyl
                    or
                    # Or if it's being transformed into something else
                    checker.check_fg("Carboxylic acid", product)
                    or checker.check_fg("Primary alcohol", product)
                    or checker.check_fg("Primary amine", product)
                    or checker.check_fg("Imine", product)
                    or checker.check_fg("Oxime", product)
                    or checker.check_fg("Acetal/Ketal", product)
                    or
                    # Or if it's part of a new ring system
                    any(
                        checker.check_ring(ring, product)
                        for ring in [
                            "benzothiazole",
                            "benzoxazole",
                            "benzimidazole",
                            "oxadiazole",
                            "thiazole",
                            "oxazole",
                            "imidazole",
                        ]
                    )
                )

                if is_formyl_manipulation or formyl_used_in_reaction:
                    found_formyl_manipulation = True
                    formyl_manipulation_reactions.append((rsmi, depth))
                    print(f"Found formyl group manipulation at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Print summary
    if found_formyl_introduction:
        print(f"Found {len(formyl_introduction_reactions)} formyl introduction reactions")
    if found_formyl_manipulation:
        print(f"Found {len(formyl_manipulation_reactions)} formyl manipulation reactions")

    # Check if we have both introduction and manipulation
    result = found_formyl_introduction and found_formyl_manipulation
    print(f"Final result: {result}")

    return result
