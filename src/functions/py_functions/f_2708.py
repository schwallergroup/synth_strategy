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
    Detects if the synthesis route uses a late-stage amide coupling as the final step.
    This looks for a reaction at depth 0 that forms an amide bond.
    """
    result = False
    min_depth_amide_coupling = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal result, min_depth_amide_coupling

        print(f"Traversing node type: {node['type']}, depth: {depth}")

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for all amide coupling reaction types
                amide_coupling_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Ester with ammonia to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Acyl chloride with ammonia to amide",
                    "Schotten-Baumann to ester",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                ]

                is_amide_coupling = False
                for reaction_type in amide_coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Detected amide coupling reaction at depth {depth}: {reaction_type}"
                        )
                        is_amide_coupling = True
                        if depth < min_depth_amide_coupling:
                            min_depth_amide_coupling = depth
                        break

                # If not identified by reaction type, check by functional groups
                if not is_amide_coupling:
                    # Alternative check: look for carboxylic acid and amine reactants
                    # and formation of a new amide bond
                    has_carboxylic_acid = False
                    has_amine = False
                    has_acyl_halide = False
                    has_ester = False

                    for reactant in reactants_smiles:
                        if checker.check_fg("Carboxylic acid", reactant):
                            has_carboxylic_acid = True
                            print(f"Found carboxylic acid in reactant: {reactant}")

                        if checker.check_fg(
                            "Primary amine", reactant
                        ) or checker.check_fg("Secondary amine", reactant):
                            has_amine = True
                            print(f"Found amine in reactant: {reactant}")

                        if checker.check_fg("Acyl halide", reactant):
                            has_acyl_halide = True
                            print(f"Found acyl halide in reactant: {reactant}")

                        if checker.check_fg("Ester", reactant):
                            has_ester = True
                            print(f"Found ester in reactant: {reactant}")

                    # Count amide groups in reactants and product
                    reactant_amide_count = 0
                    for reactant in reactants_smiles:
                        if checker.check_fg("Primary amide", reactant):
                            reactant_amide_count += 1
                            print(f"Found primary amide in reactant: {reactant}")
                        if checker.check_fg("Secondary amide", reactant):
                            reactant_amide_count += 1
                            print(f"Found secondary amide in reactant: {reactant}")
                        if checker.check_fg("Tertiary amide", reactant):
                            reactant_amide_count += 1
                            print(f"Found tertiary amide in reactant: {reactant}")

                    product_amide_count = 0
                    if checker.check_fg("Primary amide", product_smiles):
                        product_amide_count += 1
                        print(f"Found primary amide in product: {product_smiles}")
                    if checker.check_fg("Secondary amide", product_smiles):
                        product_amide_count += 1
                        print(f"Found secondary amide in product: {product_smiles}")
                    if checker.check_fg("Tertiary amide", product_smiles):
                        product_amide_count += 1
                        print(f"Found tertiary amide in product: {product_smiles}")

                    print(
                        f"Product amide count: {product_amide_count}, Reactant amide count: {reactant_amide_count}"
                    )

                    # Check if a new amide bond was formed
                    if product_amide_count > reactant_amide_count:
                        if (
                            (has_carboxylic_acid and has_amine)
                            or (has_acyl_halide and has_amine)
                            or (has_ester and has_amine)
                        ):
                            print(
                                f"Detected amide coupling at depth {depth} based on functional group analysis"
                            )
                            is_amide_coupling = True
                            if depth < min_depth_amide_coupling:
                                min_depth_amide_coupling = depth

                # If this is the final reaction (depth 0) and it's an amide coupling
                if depth == 0 and is_amide_coupling:
                    result = True

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    print("Starting traversal")
    dfs_traverse(route)
    print(f"Min depth amide coupling: {min_depth_amide_coupling}")

    # If we didn't find an amide coupling at depth 0 but found one at depth 1,
    # we'll still consider it as late-stage
    if not result and min_depth_amide_coupling <= 1:
        print(
            f"Found amide coupling at depth {min_depth_amide_coupling}, considering it late-stage"
        )
        result = True

    print(f"Final result: {result}")
    return result
