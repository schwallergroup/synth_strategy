#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    This function detects a strategy involving amide hydrolysis,
    breaking a C-N bond in an amide group.

    In retrosynthesis, we're looking for reactions where:
    - The "product" (target molecule) contains an amide
    - The "reactants" (starting materials) contain an amine and a carboxylic acid
    """
    has_amide_hydrolysis = False

    def dfs_traverse(node, depth=0):
        nonlocal has_amide_hydrolysis

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # In retrosynthesis, the "reactants" in rsmi are what we're trying to make
                # and the "products" are the starting materials
                reactants_smiles = rsmi.split(">")[0].split(".")
                products_smiles = rsmi.split(">")[-1].split(".")

                # Check if this is a hydrolysis reaction (forward direction)
                is_hydrolysis = checker.check_reaction(
                    "Hydrolysis of amides/imides/carbamates", rsmi
                ) or checker.check_reaction("Hydrogenolysis of amides/imides/carbamates", rsmi)

                if is_hydrolysis:
                    print(f"Detected hydrolysis reaction type at depth {depth}")

                    # In retrosynthesis, the product (target) should have the amide
                    # and the reactants (starting materials) should have amine and carboxylic acid
                    products_with_amide = []
                    for product in products_smiles:
                        if (
                            checker.check_fg("Primary amide", product)
                            or checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Tertiary amide", product)
                        ):
                            products_with_amide.append(product)
                            print(f"Found amide in product: {product}")

                    reactants_with_amine = []
                    reactants_with_carboxylic = []

                    for reactant in reactants_smiles:
                        if (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Tertiary amine", reactant)
                        ):
                            reactants_with_amine.append(reactant)
                            print(f"Found amine in reactant: {reactant}")

                        if checker.check_fg("Carboxylic acid", reactant):
                            reactants_with_carboxylic.append(reactant)
                            print(f"Found carboxylic acid in reactant: {reactant}")

                    # If amide is in product and amine+carboxylic acid are in reactants
                    if products_with_amide and reactants_with_amine and reactants_with_carboxylic:
                        print(
                            f"Confirmed amide hydrolysis at depth {depth}: Amide in product, amine and carboxylic acid in reactants"
                        )
                        has_amide_hydrolysis = True
                else:
                    # Check for implicit amide hydrolysis by examining functional group changes
                    print(f"Checking for implicit amide hydrolysis at depth {depth}")

                    # In retrosynthesis, check if product has amide
                    products_with_amide = [
                        p
                        for p in products_smiles
                        if checker.check_fg("Primary amide", p)
                        or checker.check_fg("Secondary amide", p)
                        or checker.check_fg("Tertiary amide", p)
                    ]

                    if products_with_amide:
                        print(f"Found amide in product for potential implicit hydrolysis")

                        # Check if reactants have both amine and carboxylic acid
                        has_amine = any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            or checker.check_fg("Tertiary amine", r)
                            for r in reactants_smiles
                        )

                        has_carboxylic = any(
                            checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
                        )

                        if has_amine and has_carboxylic:
                            print(f"Detected implicit amide hydrolysis at depth {depth}: {rsmi}")
                            print(f"Amide in product, amine and carboxylic acid in reactants")
                            has_amide_hydrolysis = True

                # Additional check for reactions that might be labeled differently
                if not has_amide_hydrolysis:
                    # Look for any reaction where amide is converted to amine and carboxylic acid
                    products_with_amide = [
                        p
                        for p in products_smiles
                        if checker.check_fg("Primary amide", p)
                        or checker.check_fg("Secondary amide", p)
                        or checker.check_fg("Tertiary amide", p)
                    ]

                    if products_with_amide:
                        has_amine = any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            or checker.check_fg("Tertiary amine", r)
                            for r in reactants_smiles
                        )

                        has_carboxylic = any(
                            checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
                        )

                        if has_amine and has_carboxylic:
                            print(
                                f"Detected functional group pattern for amide hydrolysis at depth {depth}"
                            )
                            has_amide_hydrolysis = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: has_amide_hydrolysis = {has_amide_hydrolysis}")

    return has_amide_hydrolysis
