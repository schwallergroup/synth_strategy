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
    Detects if the final step in the synthesis is an esterification reaction.
    This is a common diversification strategy.
    """
    final_step_is_esterification = False

    # Find the target molecule (root of the tree)
    target_mol = None
    if route["type"] == "mol":
        target_mol = route["smiles"]

    def dfs_traverse(node, depth=0, path_to_root=None):
        nonlocal final_step_is_esterification, target_mol

        if path_to_root is None:
            path_to_root = []

        current_path = path_to_root.copy()
        current_path.append(node)

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        # If this is a reaction node that produces the target molecule or is one step away from it
        if node["type"] == "reaction" and depth <= 1:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining potential late-stage reaction: {rsmi}")

                # Check for various esterification reactions
                if (
                    checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    or checker.check_reaction("Schotten-Baumann to ester", rsmi)
                    or checker.check_reaction(
                        "O-alkylation of carboxylic acids with diazo compounds", rsmi
                    )
                    or checker.check_reaction("Transesterification", rsmi)
                ):
                    print(f"Detected esterification reaction directly: {rsmi}")
                    final_step_is_esterification = True
                    return

                # Fallback to functional group analysis
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carboxylic acid in reactants
                acid_found = False
                alcohol_found = False

                for reactant in reactants:
                    if not reactant:
                        continue

                    print(f"Checking reactant: {reactant}")

                    if checker.check_fg("Carboxylic acid", reactant):
                        acid_found = True
                        print("Found carboxylic acid in reactant")

                    # Check for various types of alcohols
                    if (
                        checker.check_fg("Primary alcohol", reactant)
                        or checker.check_fg("Secondary alcohol", reactant)
                        or checker.check_fg("Tertiary alcohol", reactant)
                        or checker.check_fg("Aromatic alcohol", reactant)
                        or checker.check_fg("Phenol", reactant)
                    ):
                        alcohol_found = True
                        print("Found alcohol in reactant")

                # Check for ester in product
                ester_found = False
                if product:
                    print(f"Checking product: {product}")
                    if checker.check_fg("Ester", product):
                        ester_found = True
                        print("Found ester in product")

                if acid_found and alcohol_found and ester_found:
                    print("Detected esterification based on functional groups")
                    final_step_is_esterification = True
                    return

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path)

    dfs_traverse(route)
    print(f"Final result: {final_step_is_esterification}")
    return final_step_is_esterification
