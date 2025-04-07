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
    Detects if the synthesis includes a late-stage esterification
    (conversion of carboxylic acid to ester in the final steps).
    """
    esterification_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal esterification_depth

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for various esterification reactions
                esterification_reactions = [
                    "Esterification of Carboxylic Acids",
                    "O-alkylation of carboxylic acids with diazo compounds",
                    "Transesterification",
                    "Oxidative esterification of primary alcohols",
                    "Schotten-Baumann to ester",
                ]

                for reaction_type in esterification_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"{reaction_type} detected at depth {depth}, rsmi: {rsmi}")
                        if esterification_depth is None or depth < esterification_depth:
                            esterification_depth = depth
                        break

                # If no specific reaction detected, use fallback checks
                if esterification_depth is None or depth < esterification_depth:
                    # Check for carboxylic acid to ester conversion
                    has_acid_reactant = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    )
                    has_alcohol_reactant = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        for r in reactants
                    )
                    has_ester_product = checker.check_fg("Ester", product)

                    # Check for direct esterification (acid + alcohol → ester)
                    if has_acid_reactant and has_alcohol_reactant and has_ester_product:
                        print(
                            f"Carboxylic acid + alcohol to ester conversion detected at depth {depth}"
                        )
                        print(f"Reactants: {reactants}")
                        print(f"Product: {product}")
                        esterification_depth = depth

                    # Check for transesterification (ester → different ester)
                    elif any(checker.check_fg("Ester", r) for r in reactants) and has_ester_product:
                        # Make sure it's not the same ester
                        reactant_esters = [r for r in reactants if checker.check_fg("Ester", r)]
                        if any(r != product for r in reactant_esters):
                            print(f"Transesterification detected at depth {depth}")
                            esterification_depth = depth

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if esterification was found in the late stage (depth ≤ 3)
    result = esterification_depth is not None and esterification_depth <= 3
    print(f"Late stage esterification result: {result} (depth: {esterification_depth})")
    return result
