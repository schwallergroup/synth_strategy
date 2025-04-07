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
    Detects if the route involves conversion of a carboxylic acid to Weinreb amide
    followed by Grignard addition to form a ketone.
    """
    weinreb_formation = False
    weinreb_to_ketone = False

    # Track the molecules that are Weinreb amides for better connection tracking
    weinreb_amide_products = set()

    def dfs_traverse(node, depth=0):
        nonlocal weinreb_formation, weinreb_to_ketone

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for Weinreb amide to ketone conversion
                if checker.check_reaction("Ketone from Weinreb amide", rsmi):
                    weinreb_to_ketone = True
                    print(f"Found Weinreb amide to ketone conversion at depth {depth}")
                    print(f"  Reaction: {rsmi}")
                else:
                    # Manual check for Weinreb amide to ketone conversion
                    has_weinreb = False
                    has_organometallic = False
                    has_ketone = checker.check_fg("Ketone", product)

                    for reactant in reactants:
                        # Check for Weinreb amide in reactants - broader pattern
                        if (
                            checker.check_fg("Tertiary amide", reactant)
                            and "N" in reactant
                            and "O" in reactant
                        ):
                            if "CON" in reactant or "NOC" in reactant or "ONC" in reactant:
                                has_weinreb = True
                                print(f"  Found Weinreb amide in reactants: {reactant}")

                        # Check for Grignard reagent or organometallic compound
                        if (
                            "Mg" in reactant
                            or "Li" in reactant
                            or checker.check_fg("Magnesium halide", reactant)
                        ):
                            has_organometallic = True
                            print(f"  Found organometallic reagent: {reactant}")

                    if has_weinreb and has_organometallic and has_ketone:
                        weinreb_to_ketone = True
                        print(
                            f"Found Weinreb amide to ketone conversion at depth {depth} (manual check)"
                        )

                # Check for carboxylic acid to Weinreb amide conversion
                # Direct check for the reaction type
                if checker.check_reaction("Asymmetric ketones from N,N-dimethylamides", rsmi):
                    weinreb_formation = True
                    print(f"Found Weinreb amide formation reaction at depth {depth}")
                    return

                # Check for carboxylic acid to Weinreb amide conversion
                has_carboxylic_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants)
                # Broader check for N,O-dimethylhydroxylamine or similar reagent
                has_hydroxylamine = any(
                    ("[CH3:1][O:2][NH:3][CH3:4]" in r or "NOC" in r or "CON" in r or "N(C)O" in r)
                    for r in reactants
                )

                # Check if product is a Weinreb amide
                is_weinreb_amide = checker.check_fg("Tertiary amide", product) and (
                    "N" in product
                    and "O" in product
                    and (
                        "CON" in product
                        or "NOC" in product
                        or "ONC" in product
                        or "[O:2][N:3]" in product
                        or "[N:3][O:2]" in product
                    )
                )

                if has_carboxylic_acid and has_hydroxylamine and is_weinreb_amide:
                    weinreb_formation = True
                    weinreb_amide_products.add(product)
                    print(f"Found Weinreb amide formation at depth {depth}")
                    print(f"  Product Weinreb amide: {product}")

                # Alternative check: direct conversion of carboxylic acid to Weinreb amide
                if has_carboxylic_acid and is_weinreb_amide:
                    for reactant in reactants:
                        # Look for N,O-dimethylhydroxylamine or similar reagent
                        if (
                            "N" in reactant
                            and "O" in reactant
                            and (
                                "CH3" in reactant
                                or "C" in reactant
                                or "[CH3:1][O:2][NH:3][CH3:4]" in reactant
                            )
                        ):
                            weinreb_formation = True
                            weinreb_amide_products.add(product)
                            print(
                                f"Found Weinreb amide formation (alternative check) at depth {depth}"
                            )
                            print(f"  Hydroxylamine reactant: {reactant}")
                            break

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(
        f"Final result: weinreb_formation={weinreb_formation}, weinreb_to_ketone={weinreb_to_ketone}"
    )
    return weinreb_formation and weinreb_to_ketone
