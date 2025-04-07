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
    Detects etherification via tosylate strategy in the synthetic route.
    """
    etherification_found = False

    def dfs(node, depth=0):
        nonlocal etherification_found

        if etherification_found:
            return

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Williamson ether synthesis or similar reactions
                if checker.check_reaction("Williamson Ether Synthesis", rsmi):
                    # Check if one of the reactants contains tosylate
                    for reactant in reactants:
                        if checker.check_fg("Tosylate", reactant):
                            print(f"Found etherification via tosylate: {rsmi}")
                            etherification_found = True
                            break
            except Exception as e:
                print(f"Error in etherification check: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)
    print(f"Etherification via tosylate strategy detected: {etherification_found}")
    return etherification_found
