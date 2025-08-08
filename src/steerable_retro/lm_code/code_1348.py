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
    This function detects if the final step in the synthesis is an esterification.
    """
    final_step_is_esterification = False
    reaction_depths = []  # Track depths of all reactions
    esterification_depths = []  # Track depths of esterification reactions

    def dfs_traverse(node, depth=0, reaction_depth=0):
        nonlocal final_step_is_esterification

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        # Check if this is a reaction node
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(
                f"Checking reaction at depth {depth}, reaction_depth {reaction_depth}, RSMI: {rsmi}"
            )

            # Track this reaction depth
            reaction_depths.append(reaction_depth)

            try:
                # Check for various esterification reaction types
                esterification_reactions = [
                    "Esterification of Carboxylic Acids",
                    "Schotten-Baumann to ester",
                    "O-alkylation of carboxylic acids with diazo compounds",
                    "Transesterification",
                    "Oxidative esterification of primary alcohols",
                    "Acetic anhydride and alcohol to ester",
                ]

                is_esterification = False
                for reaction_type in esterification_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Detected esterification at depth {depth}, reaction_depth {reaction_depth}: {reaction_type}"
                        )
                        is_esterification = True
                        break

                # If no specific esterification reaction was detected, check for functional group changes
                if not is_esterification:
                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]

                    # Check if reactants contain carboxylic acid and product contains ester
                    reactants = reactants_part.split(".")

                    acid_in_reactants = False
                    alcohol_in_reactants = False

                    for reactant in reactants:
                        if checker.check_fg("Carboxylic acid", reactant):
                            acid_in_reactants = True
                            print(f"Found carboxylic acid in reactant: {reactant}")
                        if (
                            checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Tertiary alcohol", reactant)
                            or checker.check_fg("Aromatic alcohol", reactant)
                            or checker.check_fg("Phenol", reactant)
                        ):
                            alcohol_in_reactants = True
                            print(f"Found alcohol in reactant: {reactant}")

                    # Check for ester in product
                    ester_in_product = checker.check_fg("Ester", product_part)
                    if ester_in_product:
                        print(f"Found ester in product: {product_part}")

                    # In retrosynthesis, we might see the reverse: ester in reactants, acid in products
                    ester_in_reactants = any(checker.check_fg("Ester", r) for r in reactants)
                    acid_in_product = checker.check_fg("Carboxylic acid", product_part)

                    # Check both forward and reverse esterification patterns
                    if (acid_in_reactants and alcohol_in_reactants and ester_in_product) or (
                        ester_in_reactants and acid_in_product
                    ):
                        print(
                            f"Detected esterification at depth {depth}, reaction_depth {reaction_depth} (functional group analysis)"
                        )
                        is_esterification = True

                # If this is an esterification, track its depth
                if is_esterification:
                    esterification_depths.append(reaction_depth)
                    print(f"Added esterification at reaction_depth {reaction_depth}")

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

            # Increment reaction depth for children
            reaction_depth += 1

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, reaction_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Sort depths to find minimum
    reaction_depths.sort()
    esterification_depths.sort()

    print(f"All reaction depths: {reaction_depths}")
    print(f"Esterification depths: {esterification_depths}")

    # Check if any esterifications were found
    if esterification_depths:
        # Check if the first reaction (minimum depth) is an esterification
        if (
            esterification_depths
            and reaction_depths
            and esterification_depths[0] == reaction_depths[0]
        ):
            final_step_is_esterification = True
            print(f"Final step is esterification (reaction_depth {esterification_depths[0]})")

    print(f"Final result: {final_step_is_esterification}")
    return final_step_is_esterification
