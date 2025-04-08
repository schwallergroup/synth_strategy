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
    This function detects a synthetic strategy involving late-stage amide formation
    as the final step of the synthesis.
    """
    late_stage_amide_formation = False
    amide_formation_depths = []
    root_is_molecule = route["type"] == "mol"

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_formation

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Analyzing reaction SMILES at depth {depth}: {rsmi}")

                # Handle both > and >> formats in reaction SMILES
                if ">>" in rsmi:
                    reactants_part = rsmi.split(">>")[0]
                    product_part = rsmi.split(">>")[1]
                else:
                    parts = rsmi.split(">")
                    reactants_part = parts[0]
                    product_part = parts[-1]

                reactants_smiles = reactants_part.split(".")
                product_smiles = product_part

                # Check for amide formation reaction types using the checker
                amide_formation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with ammonia to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with primary amine to imide",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with ammonia to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                    "Carboxylic acid to amide conversion",
                    "Schotten-Baumann to ester",
                    "Acylation of secondary amines",
                    "Acylation of primary amines",
                ]

                is_amide_formation = False

                # Check if any of the amide formation reactions match
                for reaction_type in amide_formation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Detected amide formation reaction: {reaction_type} at depth {depth}"
                        )
                        is_amide_formation = True
                        break

                # If no specific reaction type matched, check for functional group changes
                if not is_amide_formation:
                    # Check for carboxylic acid, amine, and amide functional groups
                    has_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants_smiles)
                    has_primary_amine = any(
                        checker.check_fg("Primary amine", r) for r in reactants_smiles
                    )
                    has_secondary_amine = any(
                        checker.check_fg("Secondary amine", r) for r in reactants_smiles
                    )
                    has_acyl_halide = any(
                        checker.check_fg("Acyl halide", r) for r in reactants_smiles
                    )
                    has_ester = any(checker.check_fg("Ester", r) for r in reactants_smiles)
                    has_anhydride = any(checker.check_fg("Anhydride", r) for r in reactants_smiles)

                    # Check if product has amide group
                    has_primary_amide = checker.check_fg("Primary amide", product_smiles)
                    has_secondary_amide = checker.check_fg("Secondary amide", product_smiles)
                    has_tertiary_amide = checker.check_fg("Tertiary amide", product_smiles)

                    # Various amide formation pathways
                    if (
                        (has_acid or has_acyl_halide or has_ester or has_anhydride)
                        and (has_primary_amine or has_secondary_amine)
                    ) and (has_primary_amide or has_secondary_amide or has_tertiary_amide):
                        print(
                            f"Detected amide formation via functional group analysis at depth {depth}"
                        )
                        is_amide_formation = True

                if is_amide_formation:
                    amide_formation_depths.append(depth)

                    # Consider depth 0 or 1 (if root is molecule) as late-stage
                    if (depth == 0) or (root_is_molecule and depth == 1):
                        late_stage_amide_formation = True
                        print(f"Late-stage amide formation detected at depth {depth}")

            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Amide formation detected at depths: {amide_formation_depths}")
    print(f"Late-stage amide formation: {late_stage_amide_formation}")

    return late_stage_amide_formation
