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
    Detects if the final step in the synthesis is an amide bond formation.
    """
    found_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_formation

        print(f"Traversing node type: {node['type']} at depth: {depth}")

        # Check if this is a reaction node at depth 1 (final synthetic step)
        if node["type"] == "reaction" and depth == 1:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Analyzing final reaction: {rsmi}")

            # Check for amide formation using reaction checkers
            amide_formation_reactions = [
                "Carboxylic acid with primary amine to amide",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Acyl chloride with secondary amine to amide",
                "Ester with primary amine to amide",
                "Ester with secondary amine to amide",
                "Acyl chloride with ammonia to amide",
                "Ester with ammonia to amide",
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                "Schotten-Baumann_amide",
                "Acylation of primary amines",
                "Acylation of secondary amines",
            ]

            for reaction_type in amide_formation_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Detected amide formation reaction: {reaction_type}")
                    found_amide_formation = True
                    return

            # If reaction type check failed, check for functional groups
            # Check reactants for amide-forming functional groups
            has_carboxylic_acid = any(
                checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
            )
            has_acyl_halide = any(checker.check_fg("Acyl halide", r) for r in reactants_smiles)
            has_ester = any(checker.check_fg("Ester", r) for r in reactants_smiles)
            has_anhydride = any(checker.check_fg("Anhydride", r) for r in reactants_smiles)

            has_primary_amine = any(checker.check_fg("Primary amine", r) for r in reactants_smiles)
            has_secondary_amine = any(
                checker.check_fg("Secondary amine", r) for r in reactants_smiles
            )
            has_ammonia = any(r.strip() == "N" for r in reactants_smiles)

            # Check product for amide group
            has_amide_product = (
                checker.check_fg("Primary amide", product_smiles)
                or checker.check_fg("Secondary amide", product_smiles)
                or checker.check_fg("Tertiary amide", product_smiles)
            )

            print(
                f"Reactant analysis: Carboxylic acid: {has_carboxylic_acid}, Acyl halide: {has_acyl_halide}, "
                f"Ester: {has_ester}, Anhydride: {has_anhydride}"
            )
            print(
                f"Amine analysis: Primary: {has_primary_amine}, Secondary: {has_secondary_amine}, Ammonia: {has_ammonia}"
            )
            print(f"Product has amide: {has_amide_product}")

            # Check if we have both an amide-forming reagent and an amine, and the product has an amide
            if has_amide_product and (
                (has_carboxylic_acid or has_acyl_halide or has_ester or has_anhydride)
                and (has_primary_amine or has_secondary_amine or has_ammonia)
            ):
                print("Found amide formation in final step based on functional group analysis")
                found_amide_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Check if route is valid
    if not route or "type" not in route:
        print("Invalid route structure")
        return False

    # Start traversal
    dfs_traverse(route)

    return found_amide_formation
