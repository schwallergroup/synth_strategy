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
    This function detects if the synthesis includes amide formation in the late stages
    (first half of the synthesis, which corresponds to low depth in retrosynthetic tree).
    """
    amide_formation_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Extract reactants and product
                try:
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Check if product contains amide group
                    has_amide_in_product = (
                        checker.check_fg("Primary amide", product_smiles)
                        or checker.check_fg("Secondary amide", product_smiles)
                        or checker.check_fg("Tertiary amide", product_smiles)
                    )

                    # Check for amide formation using the checker function
                    amide_formation_reactions = [
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                        "Acyl chloride with ammonia to amide",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Acyl chloride with primary amine to imide",
                        "Acyl chloride with secondary amine to amide",
                        "Carboxylic acid with primary amine to amide",
                        "Ester with ammonia to amide",
                        "Ester with primary amine to amide",
                        "Ester with secondary amine to amide",
                        "Schotten-Baumann_amide",
                        "Schotten-Baumann to ester",
                        "Aminolysis of esters",
                        "Nitrile to amide",
                        "Carboxylic acid to amide conversion",
                    ]

                    is_amide_formation = False
                    for reaction_type in amide_formation_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(f"Found potential amide formation reaction: {reaction_type}")
                            print(f"Reaction SMILES: {rsmi}")

                            # Verify this is actually forming an amide by checking product
                            if has_amide_in_product:
                                print(f"Confirmed amide in product")
                                is_amide_formation = True
                                break
                            else:
                                print(f"Reaction type matches but no amide found in product")

                    # If no specific reaction type matched, check for general amide formation
                    if not is_amide_formation:
                        # Check if reactants contain necessary functional groups
                        has_acid_or_derivative = False
                        has_amine = False

                        for reactant in reactants_smiles:
                            if (
                                checker.check_fg("Carboxylic acid", reactant)
                                or checker.check_fg("Ester", reactant)
                                or checker.check_fg("Acyl halide", reactant)
                                or checker.check_fg("Nitrile", reactant)
                                or checker.check_fg("Anhydride", reactant)
                            ):
                                has_acid_or_derivative = True

                            if (
                                checker.check_fg("Primary amine", reactant)
                                or checker.check_fg("Secondary amine", reactant)
                                or checker.check_fg("Aniline", reactant)
                                or checker.check_fg("Ammonia", reactant)
                            ):
                                has_amine = True

                        # If reactants have necessary groups and product has amide, it's likely amide formation
                        if has_acid_or_derivative and has_amine and has_amide_in_product:
                            print(f"Detected amide formation from functional group analysis")
                            print(f"Reaction SMILES: {rsmi}")
                            is_amide_formation = True

                    if is_amide_formation:
                        print(f"Confirmed amide formation at depth: {depth}")
                        # If we find an amide formation, record its depth
                        if amide_formation_depth is None or depth < amide_formation_depth:
                            amide_formation_depth = depth

                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it late-stage if it's in the first half of the synthesis
    # (which means lower depth in retrosynthetic tree)
    if amide_formation_depth is not None:
        is_late_stage = amide_formation_depth < (max_depth / 2)
        print(
            f"Amide formation depth: {amide_formation_depth}, Max depth: {max_depth}, Is late stage: {is_late_stage}"
        )
        return is_late_stage

    print("No amide formation found in the synthesis route")
    return False
