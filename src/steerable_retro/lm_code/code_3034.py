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
    This function detects a linear synthesis strategy that culminates in an amide bond formation.
    """
    # Initialize tracking variables
    is_linear = True
    has_final_amide_formation = False
    min_depth_reaction = float("inf")
    final_reaction_rsmi = None

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, has_final_amide_formation, min_depth_reaction, final_reaction_rsmi

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if linear (no more than 2 reactants)
            if len(reactants_smiles) > 2:
                is_linear = False
                print(
                    f"Non-linear step detected at depth {depth} with {len(reactants_smiles)} reactants"
                )

            # Track the final reaction (minimum depth)
            if depth < min_depth_reaction:
                min_depth_reaction = depth
                final_reaction_rsmi = rsmi
                print(f"New final reaction at depth {depth}: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the final reaction is an amide formation
    if final_reaction_rsmi:
        # Check for common amide formation reactions
        amide_formation_reactions = [
            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
            "Carboxylic acid with primary amine to amide",
            "Acyl chloride with primary amine to amide (Schotten-Baumann)",
            "Acyl chloride with secondary amine to amide",
            "Ester with primary amine to amide",
            "Ester with secondary amine to amide",
            "Schotten-Baumann_amide",
        ]

        for reaction_type in amide_formation_reactions:
            if checker.check_reaction(reaction_type, final_reaction_rsmi):
                has_final_amide_formation = True
                print(f"Detected amide formation in final step: {reaction_type}")
                break

        # If no specific reaction type matched, check for amide formation by examining product
        if not has_final_amide_formation:
            product_smiles = final_reaction_rsmi.split(">")[-1]
            reactants_smiles = final_reaction_rsmi.split(">")[0].split(".")

            # Check if product has amide but reactants don't
            product_has_amide = (
                checker.check_fg("Secondary amide", product_smiles)
                or checker.check_fg("Primary amide", product_smiles)
                or checker.check_fg("Tertiary amide", product_smiles)
            )

            reactants_have_amide = False
            for reactant in reactants_smiles:
                if (
                    checker.check_fg("Secondary amide", reactant)
                    or checker.check_fg("Primary amide", reactant)
                    or checker.check_fg("Tertiary amide", reactant)
                ):
                    reactants_have_amide = True
                    break

            if product_has_amide and not reactants_have_amide:
                has_final_amide_formation = True
                print("Detected amide formation in final step (by functional group analysis)")

    # Check if the strategy is present
    strategy_present = is_linear and has_final_amide_formation
    if strategy_present:
        print("Detected linear synthesis with final amide formation")
    else:
        if not is_linear:
            print("Synthesis is not linear")
        if not has_final_amide_formation:
            print("Final step does not form an amide bond")

    return strategy_present
