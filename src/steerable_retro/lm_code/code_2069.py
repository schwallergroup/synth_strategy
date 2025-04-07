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
    This function detects if the synthesis involves early-stage pyridone ring formation.
    """
    pyridone_formed = False
    high_depth_threshold = 3  # Depth >= 3 is considered early stage

    def dfs_traverse(node, depth=0):
        nonlocal pyridone_formed

        if node["type"] == "reaction" and depth >= high_depth_threshold:
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                # Check if pyridone is in product but not in any reactants (true formation)
                product_has_pyridone = checker.check_ring("pyridone", product)
                reactants_have_pyridone = any(
                    checker.check_ring("pyridone", reactant) for reactant in reactants
                )

                # Check for pyridine to pyridone conversion
                reactants_have_pyridine = any(
                    checker.check_ring("pyridine", reactant) for reactant in reactants
                )

                if product_has_pyridone and not reactants_have_pyridone:
                    print(f"Found pyridone formation at depth {depth}")
                    pyridone_formed = True
                elif (
                    product_has_pyridone and reactants_have_pyridine and not reactants_have_pyridone
                ):
                    print(f"Found pyridine to pyridone conversion at depth {depth}")
                    pyridone_formed = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return pyridone_formed
