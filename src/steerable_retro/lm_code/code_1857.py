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
    This function detects routes with amide coupling as one of the final steps (depth 0 or 1).
    """
    has_late_stage_amide_coupling = False

    def dfs_traverse(node, current_depth=0):
        nonlocal has_late_stage_amide_coupling

        if node["type"] == "reaction":
            # Only consider reactions at depth 0 or 1 (late stage)
            if current_depth <= 1:
                # Extract reaction SMILES
                try:
                    rsmi = node["metadata"]["rsmi"]

                    # Check for amide coupling reactions directly
                    amide_coupling_reactions = [
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                        "Acyl chloride with ammonia to amide",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Acyl chloride with secondary amine to amide",
                        "Carboxylic acid with primary amine to amide",
                        "Ester with ammonia to amide",
                        "Ester with primary amine to amide",
                        "Ester with secondary amine to amide",
                        "Schotten-Baumann to ester",
                        "Schotten-Baumann_amide",
                    ]

                    for reaction_type in amide_coupling_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(
                                f"Detected late-stage amide coupling reaction: {reaction_type} at depth {current_depth}"
                            )
                            has_late_stage_amide_coupling = True
                            return

                    # If no specific reaction type matched, check for functional group changes
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Check if any reactant has amine (primary or secondary)
                    has_primary_amine = any(
                        checker.check_fg("Primary amine", r) for r in reactants_smiles
                    )
                    has_secondary_amine = any(
                        checker.check_fg("Secondary amine", r) for r in reactants_smiles
                    )
                    has_amine = has_primary_amine or has_secondary_amine

                    # Check if any reactant has carboxylic acid or acyl halide
                    has_carboxylic = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
                    )
                    has_acyl_halide = any(
                        checker.check_fg("Acyl halide", r) for r in reactants_smiles
                    )

                    # Check if product has amide
                    has_primary_amide = checker.check_fg("Primary amide", product_smiles)
                    has_secondary_amide = checker.check_fg("Secondary amide", product_smiles)
                    has_tertiary_amide = checker.check_fg("Tertiary amide", product_smiles)
                    has_amide_product = (
                        has_primary_amide or has_secondary_amide or has_tertiary_amide
                    )

                    if has_amine and (has_carboxylic or has_acyl_halide) and has_amide_product:
                        print(
                            f"Detected late-stage amide coupling at depth {current_depth} through functional group analysis"
                        )
                        has_late_stage_amide_coupling = True
                except Exception as e:
                    print(f"Error analyzing reaction: {e}")

        # Traverse children (depth increases as we go deeper in the retrosynthetic tree)
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_late_stage_amide_coupling
