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
    This function detects amide formation from acid chloride or other acylating agents
    as one of the final steps in the synthesis (within last 3 steps).
    """
    late_stage_amide_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_formation_found

        if (
            node["type"] == "reaction" and depth <= 2
        ):  # Check if it's a late-stage step (depth 0-2)
            try:
                # Extract reactants and product from reaction SMILES
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is an amide formation reaction using various reaction patterns
                is_amide_formation = False

                # Check for acyl halide-based amide formation
                if (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction("Schotten-Baumann to ester", rsmi)
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acyl chloride with secondary amine to amide", rsmi
                    )
                    or checker.check_reaction(
                        "Acyl chloride with ammonia to amide", rsmi
                    )
                    or checker.check_reaction("Acylation of primary amines", rsmi)
                    or checker.check_reaction("Acylation of secondary amines", rsmi)
                ):
                    is_amide_formation = True
                    print("Detected acyl halide-based amide formation reaction")

                # Check for carboxylic acid-based amide formation
                if not is_amide_formation:
                    if checker.check_reaction(
                        "Carboxylic acid with primary amine to amide", rsmi
                    ) or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    ):
                        is_amide_formation = True
                        print("Detected carboxylic acid-based amide formation reaction")

                # Check for ester-based amide formation
                if not is_amide_formation:
                    if (
                        checker.check_reaction("Ester with ammonia to amide", rsmi)
                        or checker.check_reaction(
                            "Ester with primary amine to amide", rsmi
                        )
                        or checker.check_reaction(
                            "Ester with secondary amine to amide", rsmi
                        )
                        or checker.check_reaction("Aminolysis of esters", rsmi)
                    ):
                        is_amide_formation = True
                        print("Detected ester-based amide formation reaction")

                # Check for acylating agent in reactants
                acylating_agent_found = False
                for reactant in reactants_smiles:
                    if checker.check_fg("Acyl halide", reactant):
                        print(f"Found acyl halide in reactant: {reactant}")
                        acylating_agent_found = True
                        break
                    elif checker.check_fg("Carboxylic acid", reactant):
                        print(f"Found carboxylic acid in reactant: {reactant}")
                        acylating_agent_found = True
                        break
                    elif checker.check_fg("Ester", reactant):
                        print(f"Found ester in reactant: {reactant}")
                        acylating_agent_found = True
                        break
                    elif checker.check_fg("Anhydride", reactant):
                        print(f"Found anhydride in reactant: {reactant}")
                        acylating_agent_found = True
                        break

                # Check for amide in product
                amide_in_product = (
                    checker.check_fg("Primary amide", product_smiles)
                    or checker.check_fg("Secondary amide", product_smiles)
                    or checker.check_fg("Tertiary amide", product_smiles)
                )

                if amide_in_product:
                    print(f"Found amide in product: {product_smiles}")

                # Verify conditions are met - either the reaction is recognized as amide formation
                # or we have both an acylating agent and an amide product
                if (is_amide_formation and amide_in_product) or (
                    acylating_agent_found and amide_in_product
                ):
                    print(f"Late-stage amide formation detected at depth {depth}")
                    late_stage_amide_formation_found = True
                else:
                    print(
                        f"Not a late-stage amide formation: acylating_agent={acylating_agent_found}, amide={amide_in_product}, amide_formation_reaction={is_amide_formation}"
                    )

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return late_stage_amide_formation_found
