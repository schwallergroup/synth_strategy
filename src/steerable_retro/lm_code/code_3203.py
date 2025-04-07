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
    This function detects a synthetic strategy involving esterification of a carboxylic acid.
    In retrosynthesis, this includes hydrolysis reactions viewed in reverse.
    """
    has_esterification = False

    def dfs_traverse(node, depth=0):
        nonlocal has_esterification

        # Print node information for debugging
        indent = "  " * depth
        if node["type"] == "mol":
            print(f"{indent}Molecule: {node['smiles']}")
        elif node["type"] == "reaction":
            print(f"{indent}Reaction node")

            # Extract reactants and product from reaction SMILES
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"{indent}Reaction SMILES: {rsmi}")

                # Check if this is an esterification reaction using the checker
                # Note: In retrosynthesis, hydrolysis reactions are viewed as esterification
                if (
                    checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    or checker.check_reaction("Transesterification", rsmi)
                    or checker.check_reaction("Schotten-Baumann to ester", rsmi)
                    or checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                    )
                ):
                    print(f"{indent}Found esterification reaction (or its reverse)!")
                    has_esterification = True
                else:
                    # Manual check for esterification pattern
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check if any reactant has carboxylic acid and product has ester (forward esterification)
                        has_acid_in_reactants = False
                        for reactant in reactants:
                            if checker.check_fg("Carboxylic acid", reactant):
                                has_acid_in_reactants = True
                                print(f"{indent}  Found carboxylic acid in reactants: {reactant}")
                                break

                        has_ester_in_reactants = False
                        for reactant in reactants:
                            if checker.check_fg("Ester", reactant):
                                has_ester_in_reactants = True
                                print(f"{indent}  Found ester in reactants: {reactant}")
                                break

                        has_acid_in_product = checker.check_fg("Carboxylic acid", product)
                        has_ester_in_product = checker.check_fg("Ester", product)

                        # Check for forward esterification (acid → ester)
                        if has_acid_in_reactants and has_ester_in_product:
                            print(f"{indent}  Found ester in product: {product}")
                            print(f"{indent}  This appears to be a forward esterification reaction")
                            has_esterification = True

                        # Check for reverse esterification (ester → acid, which is hydrolysis)
                        elif has_ester_in_reactants and has_acid_in_product:
                            print(f"{indent}  Found carboxylic acid in product: {product}")
                            print(
                                f"{indent}  This appears to be a hydrolysis reaction (reverse esterification)"
                            )
                            has_esterification = True

                        # Check for other ester-forming reactions
                        elif not has_esterification:
                            has_acyl_halide = False
                            has_alcohol = False

                            for reactant in reactants:
                                if checker.check_fg("Acyl halide", reactant):
                                    has_acyl_halide = True
                                    print(f"{indent}  Found acyl halide in reactants: {reactant}")
                                if (
                                    checker.check_fg("Primary alcohol", reactant)
                                    or checker.check_fg("Secondary alcohol", reactant)
                                    or checker.check_fg("Tertiary alcohol", reactant)
                                ):
                                    has_alcohol = True
                                    print(f"{indent}  Found alcohol in reactants: {reactant}")

                            if has_acyl_halide and has_alcohol and has_ester_in_product:
                                print(f"{indent}  This appears to be an acyl halide esterification")
                                has_esterification = True
                    except Exception as e:
                        print(f"{indent}Error analyzing reaction components: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting traversal of synthetic route")
    dfs_traverse(route)

    print(f"Esterification strategy found: {has_esterification}")
    return has_esterification
