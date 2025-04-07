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
    This function detects a synthetic strategy involving carbonyl group formation.
    """
    carbonyl_formed = False

    # List of carbonyl functional groups to check
    carbonyl_fgs = [
        "Aldehyde",
        "Ketone",
        "Carboxylic acid",
        "Ester",
        "Primary amide",
        "Secondary amide",
        "Tertiary amide",
        "Anhydride",
        "Acyl halide",
    ]

    # List of reactions that typically form carbonyl groups
    carbonyl_forming_reactions = [
        "Oxidation of aldehydes to carboxylic acids",
        "Oxidation of alcohols to carboxylic acids",
        "Oxidation of ketone to carboxylic acid",
        "Oxidation of nitrile to carboxylic acid",
        "Oxidation of alkene to carboxylic acid",
        "Oxidation of aldehydes and ketones to alcohols",
        "Oxidation of alcohol to aldehyde",
        "Oxidation of alcohol and aldehyde to ester",
        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
        "Hydration of alkyne to ketone",
        "Hydration of alkyne to aldehyde",
        "Ketone from Weinreb amide",
        "Friedel-Crafts acylation",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Acylation of secondary amines with anhydrides",
        "Acylation of secondary amines",
        "Acylation of primary amines",
        "Ester and halide to ketone",
        "Asymmetric ketones from N,N-dimethylamides",
        "Carboxylic acid from Li and CO2",
        "Ketone from Li and CO2",
        "Ketone from Li, Grignard and CO2",
        "Ketone from Li, halide and CO2",
        "Ketone from Grignard and CO2",
        "Grignard from nitrile to ketone",
        "Nef reaction (nitro to ketone)",
        "Acylation of olefines by aldehydes",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal carbonyl_formed

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a carbonyl-forming reaction type
                if any(checker.check_reaction(rxn, rsmi) for rxn in carbonyl_forming_reactions):
                    print(f"Carbonyl-forming reaction detected at depth {depth}: {rsmi}")
                    carbonyl_formed = True

                # Check for carbonyl groups in product and reactants
                product_has_carbonyl = any(checker.check_fg(fg, product) for fg in carbonyl_fgs)
                reactants_have_carbonyl = any(
                    any(checker.check_fg(fg, r) for fg in carbonyl_fgs) for r in reactants
                )

                # Count specific carbonyl groups in product
                product_carbonyl_count = sum(
                    1 for fg in carbonyl_fgs if checker.check_fg(fg, product)
                )

                # Count specific carbonyl groups in reactants
                reactant_carbonyl_count = sum(
                    sum(1 for fg in carbonyl_fgs if checker.check_fg(fg, r)) for r in reactants
                )

                # In retrosynthesis, if reactants have carbonyl but product doesn't, or reactants have more carbonyls
                # This indicates carbonyl formation in the forward direction
                if (reactants_have_carbonyl and not product_has_carbonyl) or (
                    reactant_carbonyl_count > product_carbonyl_count
                ):
                    print(f"Carbonyl formation detected at depth {depth}: {rsmi}")
                    print(
                        f"Product carbonyl count: {product_carbonyl_count}, Reactant carbonyl count: {reactant_carbonyl_count}"
                    )
                    carbonyl_formed = True

                # Check for specific transformations that create carbonyls
                if not carbonyl_formed:
                    # Check for alcohol oxidation to aldehyde/ketone (in retrosynthesis)
                    reactants_have_aldehyde_or_ketone = any(
                        checker.check_fg("Aldehyde", r) or checker.check_fg("Ketone", r)
                        for r in reactants
                    )
                    product_has_alcohol = (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                    )

                    if reactants_have_aldehyde_or_ketone and product_has_alcohol:
                        print(f"Alcohol oxidation to carbonyl detected at depth {depth}: {rsmi}")
                        carbonyl_formed = True

                    # Check for alkyne hydration to aldehyde/ketone (in retrosynthesis)
                    product_has_alkyne = checker.check_fg("Alkyne", product)
                    if reactants_have_aldehyde_or_ketone and product_has_alkyne:
                        print(f"Alkyne hydration to carbonyl detected at depth {depth}: {rsmi}")
                        carbonyl_formed = True

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if carbonyl_formed:
        print("Carbonyl formation strategy detected")
    else:
        print("No carbonyl formation strategy detected")

    return carbonyl_formed
