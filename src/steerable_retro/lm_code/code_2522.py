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
    This function detects if the synthetic route involves a chain of at least 4 distinct
    functional group transformations.
    """
    # Track the chain of transformations
    transformation_chain = []

    # Define common functional group interconversions to check
    fg_transformations = [
        # Format: (reactant_fg, product_fg, reaction_type)
        ("Aldehyde", "Alcohol", "Reduction of aldehydes and ketones to alcohols"),
        ("Ketone", "Alcohol", "Reduction of aldehydes and ketones to alcohols"),
        ("Carboxylic acid", "Ester", "Esterification of Carboxylic Acids"),
        (
            "Ester",
            "Carboxylic acid",
            "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
        ),
        ("Ester", "Carboxylic acid", "Ester saponification (methyl deprotection)"),
        ("Ester", "Carboxylic acid", "Ester saponification (alkyl deprotection)"),
        ("Alcohol", "Halide", "Alcohol to chloride_SOCl2"),
        ("Alcohol", "Halide", "Alcohol to bromide"),
        ("Alcohol", "Halide", "Appel reaction"),
        ("Alcohol", "Aldehyde", "Oxidation of Alcohols to Aldehydes and Ketones"),
        ("Alcohol", "Ketone", "Oxidation of Alcohols to Aldehydes and Ketones"),
        ("Alkene", "Alcohol", "Alkene to diol"),
        ("Alkene", "Alcohol", "anti-Markovnikov alkene hydration to alcohol"),
        ("Alkene", "Alcohol", "Markovnikov alkene hydration to alcohol"),
        ("Nitrile", "Amine", "Reduction of nitrile to amine"),
        ("Nitrile", "Amide", "Nitrile to amide"),
        ("Nitro group", "Amine", "Reduction of nitro groups to amines"),
        ("Carboxylic acid", "Amide", "Carboxylic acid to amide conversion"),
        ("Carboxylic acid", "Acyl halide", "Acyl chloride with ammonia to amide"),
        ("Acyl halide", "Amide", "Acyl chloride with primary amine to amide (Schotten-Baumann)"),
        ("Acyl halide", "Ester", "Schotten-Baumann to ester"),
        ("Aldehyde", "Oxime", "Ketone/aldehyde to hydrazone"),
        ("Ketone", "Oxime", "Ketone/aldehyde to hydrazone"),
        ("Alcohol", "Ether", "Williamson Ether Synthesis"),
        ("Alcohol", "Ester", "Acetic anhydride and alcohol to ester"),
        ("Alcohol", "Halide", "PBr3 and alcohol to alkyl bromide"),
        ("Alcohol", "Triflate", "Alcohol to triflate conversion"),
        ("Amine", "Azide", "Amine to azide"),
        ("Azide", "Amine", "Azide to amine reduction (Staudinger)"),
        ("Boronic acid", "Boronic ester", "Preparation of boronic ethers"),
        ("Boronic ester", "Boronic acid", "Preparation of boronic acids from boronic ether"),
        ("Carboxylic acid", "Alcohol", "Reduction of carboxylic acid to primary alcohol"),
        ("Alcohol", "Carboxylic acid", "Oxidation of alcohol to carboxylic acid"),
        ("Alcohol", "Alkyl halide", "Primary alkyl halide to alcohol"),
        ("Aldehyde", "Carboxylic acid", "Oxidation of aldehydes to carboxylic acids"),
        (
            "Primary alcohol",
            "Aldehyde",
            "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
        ),
        (
            "Secondary alcohol",
            "Ketone",
            "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
        ),
        ("Primary amine", "Secondary amine", "N-alkylation of primary amines with alkyl halides"),
        (
            "Secondary amine",
            "Tertiary amine",
            "N-alkylation of secondary amines with alkyl halides",
        ),
        ("Primary halide", "Primary amine", "N-alkylation of primary amines with alkyl halides"),
        (
            "Primary halide",
            "Secondary amine",
            "N-alkylation of secondary amines with alkyl halides",
        ),
    ]

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Analyzing reaction: {rsmi}")

            # Extract reactants and products
            reactants_part = rsmi.split(">")[0]
            agents_part = rsmi.split(">")[1]
            products_part = rsmi.split(">")[2]

            reactants = reactants_part.split(".")
            products = products_part.split(".")

            # Check for each functional group transformation
            for reactant_fg, product_fg, reaction_type in fg_transformations:
                # Check if this reaction matches the expected type
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Reaction matches type: {reaction_type}")

                    # Check if reactants have the expected functional group
                    reactant_has_fg = False
                    for reactant in reactants:
                        if checker.check_fg(reactant_fg, reactant):
                            reactant_has_fg = True
                            print(f"Found {reactant_fg} in reactant: {reactant}")
                            break

                    # Check if products have the expected functional group
                    product_has_fg = False
                    for product in products:
                        if checker.check_fg(product_fg, product):
                            product_has_fg = True
                            print(f"Found {product_fg} in product: {product}")
                            break

                    if reactant_has_fg and product_has_fg:
                        transformation = (reactant_fg, product_fg, reaction_type)
                        transformation_chain.append(transformation)
                        print(
                            f"Verified transformation: {reactant_fg} -> {product_fg} via {reaction_type}"
                        )
                        break  # Found a valid transformation for this reaction

                # Special case for aldehyde to oxime transformation
                elif reactant_fg == "Aldehyde" and product_fg == "Oxime":
                    # Check if reactants contain aldehyde and products contain oxime
                    reactant_has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants)
                    product_has_oxime = any(checker.check_fg("Oxime", p) for p in products)

                    if reactant_has_aldehyde and product_has_oxime:
                        print(f"Found special case: Aldehyde -> Oxime")
                        transformation = (reactant_fg, product_fg, "Aldehyde to Oxime")
                        transformation_chain.append(transformation)
                        break

                # Special case for alcohol to aldehyde oxidation
                elif reactant_fg == "Alcohol" and product_fg == "Aldehyde":
                    # Check if reactants contain alcohol and products contain aldehyde
                    reactant_has_alcohol = any(checker.check_fg("Alcohol", r) for r in reactants)
                    product_has_aldehyde = any(checker.check_fg("Aldehyde", p) for p in products)

                    if reactant_has_alcohol and product_has_aldehyde:
                        print(f"Found special case: Alcohol -> Aldehyde")
                        transformation = (reactant_fg, product_fg, "Alcohol to Aldehyde")
                        transformation_chain.append(transformation)
                        break

                # Special case for primary amine alkylation
                elif reactant_fg == "Primary amine" and product_fg == "Secondary amine":
                    # Check if reactants contain primary amine and products contain secondary amine
                    reactant_has_primary_amine = any(
                        checker.check_fg("Primary amine", r) for r in reactants
                    )
                    product_has_secondary_amine = any(
                        checker.check_fg("Secondary amine", p) for p in products
                    )

                    if reactant_has_primary_amine and product_has_secondary_amine:
                        print(f"Found special case: Primary amine -> Secondary amine")
                        transformation = (reactant_fg, product_fg, "Primary amine alkylation")
                        transformation_chain.append(transformation)
                        break

                # Special case for carboxylic acid to alcohol reduction
                elif reactant_fg == "Carboxylic acid" and product_fg == "Alcohol":
                    # Check if reactants contain carboxylic acid and products contain alcohol
                    reactant_has_carboxylic_acid = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    )
                    product_has_alcohol = any(checker.check_fg("Alcohol", p) for p in products)

                    if reactant_has_carboxylic_acid and product_has_alcohol:
                        print(f"Found special case: Carboxylic acid -> Alcohol")
                        transformation = (reactant_fg, product_fg, "Carboxylic acid reduction")
                        transformation_chain.append(transformation)
                        break

        # Traverse children (retrosynthetic direction)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least 4 distinct transformations
    distinct_transformations = set()
    for reactant_fg, product_fg, _ in transformation_chain:
        distinct_transformations.add((reactant_fg, product_fg))

    # Print all transformations found
    print(f"All transformations found: {transformation_chain}")

    result = len(distinct_transformations) >= 4
    print(
        f"Functional group interconversion chain detected: {result} (count: {len(distinct_transformations)})"
    )
    print(f"Distinct transformations: {distinct_transformations}")
    return result
