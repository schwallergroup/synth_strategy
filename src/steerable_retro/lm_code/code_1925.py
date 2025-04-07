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
    This function detects a late-stage amide coupling strategy where a carboxylic acid
    and an amine/sulfonamide are connected in the final step of the synthesis.
    """
    is_late_stage_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal is_late_stage_amide_coupling

        print(f"Traversing node type: {node['type']} at depth: {depth}")

        if node["type"] == "reaction" and depth <= 1:  # Final or penultimate reaction (late-stage)
            if "metadata" in node and "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Analyzing reaction at depth {depth}: {rsmi}")

                    # Check if this is an amide coupling reaction
                    amide_coupling_reactions = [
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                        "Carboxylic acid with primary amine to amide",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Acyl chloride with secondary amine to amide",
                        "Schotten-Baumann_amide",
                    ]

                    is_amide_coupling = False
                    for reaction_type in amide_coupling_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(f"Detected amide coupling reaction: {reaction_type}")
                            is_amide_coupling = True
                            break

                    # If not detected by reaction type, check for functional groups
                    if not is_amide_coupling:
                        print("Checking functional groups for amide coupling...")

                        # Check if reactants contain carboxylic acid and amine/sulfonamide
                        has_acid = False
                        has_acyl_halide = False
                        has_amine = False
                        has_sulfonamide = False

                        for reactant in reactants:
                            if checker.check_fg("Carboxylic acid", reactant):
                                print(f"Found carboxylic acid in reactant: {reactant}")
                                has_acid = True
                            if checker.check_fg("Acyl halide", reactant):
                                print(f"Found acyl halide in reactant: {reactant}")
                                has_acyl_halide = True
                            if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                                "Secondary amine", reactant
                            ):
                                print(f"Found amine in reactant: {reactant}")
                                has_amine = True
                            if checker.check_fg("Sulfonamide", reactant):
                                print(f"Found sulfonamide in reactant: {reactant}")
                                has_sulfonamide = True

                        # Check if product has amide bond
                        has_amide = (
                            checker.check_fg("Primary amide", product)
                            or checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Tertiary amide", product)
                        )

                        if has_amide:
                            print(f"Found amide in product: {product}")
                            if (has_acid or has_acyl_halide) and (has_amine or has_sulfonamide):
                                print("Functional group analysis confirms amide coupling")
                                is_amide_coupling = True

                    if is_amide_coupling and depth <= 1:
                        is_late_stage_amide_coupling = True
                        print(f"Confirmed late-stage amide coupling at depth {depth}")
                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    print(f"Final result: {is_late_stage_amide_coupling}")
    return is_late_stage_amide_coupling
