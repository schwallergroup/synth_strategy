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
    This function detects if the synthesis involves a trifluoromethyl-containing
    heterocycle as a key building block.
    """
    found_cf3_heterocycle = False

    # List of common heterocycles to check
    heterocycles = [
        "pyridine",
        "pyrrole",
        "furan",
        "thiophene",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
    ]

    def dfs_traverse(node):
        nonlocal found_cf3_heterocycle

        if node["type"] == "mol":
            smiles = node["smiles"]
            try:
                # Check if molecule contains CF3 group
                has_cf3 = checker.check_fg("Trifluoro group", smiles)

                if has_cf3:
                    print(f"Found molecule with CF3 group: {smiles}")

                    # Check if molecule also contains a heterocycle
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, smiles):
                            print(
                                f"Found trifluoromethyl-containing heterocycle ({heterocycle}): {smiles}"
                            )
                            found_cf3_heterocycle = True
                            break
            except Exception as e:
                print(f"Error checking molecule {smiles}: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_cf3_heterocycle
