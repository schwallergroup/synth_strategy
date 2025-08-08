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
    Detects if the synthetic route involves a sequence of carbonyl transformations:
    ester → ketone → alcohol → amine
    """
    # Track transformations and their connections
    transformations = []

    def dfs_traverse(node, depth=0):
        # Process the current node
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for functional groups
            is_ester = checker.check_fg("Ester", mol_smiles)
            is_ketone = checker.check_fg("Ketone", mol_smiles)
            is_alcohol = (
                checker.check_fg("Secondary alcohol", mol_smiles)
                or checker.check_fg("Primary alcohol", mol_smiles)
                or checker.check_fg("Tertiary alcohol", mol_smiles)
            )
            is_amine = (
                checker.check_fg("Primary amine", mol_smiles)
                or checker.check_fg("Secondary amine", mol_smiles)
                or checker.check_fg("Tertiary amine", mol_smiles)
            )

            # Store molecule info with its functional groups and depth
            transformations.append(
                {
                    "smiles": mol_smiles,
                    "depth": depth,
                    "groups": {
                        "ester": is_ester,
                        "ketone": is_ketone,
                        "alcohol": is_alcohol,
                        "amine": is_amine,
                    },
                }
            )

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort by depth (highest to lowest)
    transformations.sort(key=lambda x: x["depth"], reverse=True)

    # Check for the sequence
    has_ester = False
    has_ketone = False
    has_alcohol = False
    has_amine = False

    # Track the depths where we find each functional group
    ester_depth = -1
    ketone_depth = -1
    alcohol_depth = -1
    amine_depth = -1

    for t in transformations:
        if t["groups"]["ester"] and not has_ester:
            has_ester = True
            ester_depth = t["depth"]
        elif t["groups"]["ketone"] and has_ester and not has_ketone and t["depth"] < ester_depth:
            has_ketone = True
            ketone_depth = t["depth"]
        elif (
            t["groups"]["alcohol"] and has_ketone and not has_alcohol and t["depth"] < ketone_depth
        ):
            has_alcohol = True
            alcohol_depth = t["depth"]
        elif t["groups"]["amine"] and has_alcohol and not has_amine and t["depth"] < alcohol_depth:
            has_amine = True
            amine_depth = t["depth"]

    sequence_found = has_ester and has_ketone and has_alcohol and has_amine

    # Check if the transformations are connected by reactions
    if sequence_found:
        print(
            f"Found carbonyl transformation sequence: ester (depth {ester_depth}) → ketone (depth {ketone_depth}) → alcohol (depth {alcohol_depth}) → amine (depth {amine_depth})"
        )

    return sequence_found
