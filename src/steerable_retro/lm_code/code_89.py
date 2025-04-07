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
    This function detects if the synthesis route employs a late-stage amide formation strategy,
    specifically looking for carboxylic acid to amide transformation in the final steps.
    """
    amide_formation_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an amide formation reaction
                is_amide_formation = False

                # First check if this is a known amide formation reaction type
                amide_reaction_types = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Carboxylic acid with primary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                ]

                for reaction_type in amide_reaction_types:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found amide formation reaction: {reaction_type} at depth {depth}")
                        is_amide_formation = True
                        break

                # If not a known reaction type, check for functional group transformation
                if not is_amide_formation:
                    # Check for carboxylic acid in reactants
                    reactant_has_acid = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    )

                    # Check for amide in product
                    product_has_amide = any(
                        checker.check_fg(amide_type, product)
                        for amide_type in ["Primary amide", "Secondary amide", "Tertiary amide"]
                    )

                    if reactant_has_acid and product_has_amide:
                        print(f"Found carboxylic acid to amide transformation at depth {depth}")
                        is_amide_formation = True

                if is_amide_formation:
                    # In retrosynthesis, lower depth means later stage in forward synthesis
                    if amide_formation_depth is None or depth < amide_formation_depth:
                        amide_formation_depth = depth
                        print(f"Updated amide formation depth to {depth}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # In retrosynthesis, late stage means close to the target (low depth)
    # Define late stage as occurring in the first third of the max possible depth
    is_late_stage = amide_formation_depth is not None and amide_formation_depth <= max_depth / 3

    print(
        f"Amide formation depth: {amide_formation_depth}, Max depth: {max_depth}, Is late stage: {is_late_stage}"
    )
    return is_late_stage
