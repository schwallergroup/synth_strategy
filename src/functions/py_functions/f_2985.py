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
    This function detects if the synthesis targets a hydroxamic acid as the final product.
    """
    hydroxamic_acid_present = False

    def dfs_traverse(node, depth=0):
        nonlocal hydroxamic_acid_present

        if node["type"] == "mol":
            # Check if this is the target molecule (root of the tree)
            if depth == 0:
                print(f"Checking final product: {node['smiles']}")

                # Use RDKit to check for hydroxamic acid pattern
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Hydroxamic acid pattern: C(=O)N-OH
                    hydroxamic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3][OX2H]")
                    if mol.HasSubstructMatch(hydroxamic_pattern):
                        print(
                            f"Hydroxamic acid detected in final product using SMARTS pattern"
                        )
                        hydroxamic_acid_present = True
                    # Alternative check using SMILES substring
                    elif "C(=O)NO" in node["smiles"]:
                        print(
                            f"Hydroxamic acid detected in final product using SMILES substring"
                        )
                        hydroxamic_acid_present = True
                    # Check for combination of amide and hydroxyl
                    elif (
                        checker.check_fg("Primary amide", node["smiles"])
                        and "OH" in node["smiles"]
                    ):
                        print(f"Potential hydroxamic acid detected (amide + hydroxyl)")
                        hydroxamic_acid_present = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return hydroxamic_acid_present
