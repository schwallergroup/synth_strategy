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
    Detects if the synthesis route involves a late-stage esterification of a carboxylic acid.
    Late stage means at depth 0 or 1 (final or penultimate step).
    """
    esterification_found = False
    esterification_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal esterification_found, esterification_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}")
                print(f"Reactants: {reactants_smiles}")
                print(f"Product: {product_smiles}")

                # Check for esterification pattern (forward direction)
                reactant_has_carboxylic = any(
                    checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
                )
                product_has_ester = checker.check_fg("Ester", product_smiles)

                # Check for hydrolysis pattern (reverse direction in retrosynthesis)
                reactant_has_ester = any(checker.check_fg("Ester", r) for r in reactants_smiles)
                product_has_carboxylic = checker.check_fg("Carboxylic acid", product_smiles)

                print(f"Reactant has carboxylic acid: {reactant_has_carboxylic}")
                print(f"Product has ester: {product_has_ester}")
                print(f"Reactant has ester: {reactant_has_ester}")
                print(f"Product has carboxylic acid: {product_has_carboxylic}")

                # Check for various esterification reactions
                is_esterification = checker.check_reaction(
                    "Esterification of Carboxylic Acids", rsmi
                )
                is_transesterification = checker.check_reaction("Transesterification", rsmi)
                is_o_alkylation = checker.check_reaction(
                    "O-alkylation of carboxylic acids with diazo compounds", rsmi
                )

                # Check for hydrolysis (reverse of esterification in retrosynthesis)
                is_hydrolysis = checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                )
                is_saponification_methyl = checker.check_reaction(
                    "Ester saponification (methyl deprotection)", rsmi
                )
                is_saponification_alkyl = checker.check_reaction(
                    "Ester saponification (alkyl deprotection)", rsmi
                )

                print(f"Is esterification reaction: {is_esterification}")
                print(f"Is transesterification: {is_transesterification}")
                print(f"Is O-alkylation: {is_o_alkylation}")
                print(f"Is hydrolysis: {is_hydrolysis}")
                print(f"Is saponification (methyl): {is_saponification_methyl}")
                print(f"Is saponification (alkyl): {is_saponification_alkyl}")

                # Forward direction: carboxylic acid -> ester
                forward_esterification = (
                    reactant_has_carboxylic
                    and product_has_ester
                    and (is_esterification or is_transesterification or is_o_alkylation)
                )

                # Reverse direction: ester -> carboxylic acid (in retrosynthesis, this is esterification)
                reverse_esterification = (
                    reactant_has_ester
                    and product_has_carboxylic
                    and (is_hydrolysis or is_saponification_methyl or is_saponification_alkyl)
                )

                # If reaction checkers fail, fall back to functional group transformation
                fg_transformation = (reactant_has_carboxylic and product_has_ester) or (
                    reactant_has_ester and product_has_carboxylic
                )

                if forward_esterification:
                    print(f"Forward esterification detected at depth {depth}")
                    esterification_found = True
                    esterification_depth = depth
                elif reverse_esterification:
                    print(f"Reverse esterification (hydrolysis) detected at depth {depth}")
                    esterification_found = True
                    esterification_depth = depth
                elif fg_transformation:
                    print(
                        f"Esterification/hydrolysis detected by functional group transformation at depth {depth}"
                    )
                    esterification_found = True
                    esterification_depth = depth
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if esterification was found at depth 0 or 1 (late stage)
    result = (
        esterification_found and (esterification_depth is not None) and (esterification_depth <= 1)
    )
    print(f"Late-stage esterification strategy detected: {result}")
    return result
