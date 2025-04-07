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
    This function detects a synthetic strategy involving a pyrazole heterocycle scaffold.
    """
    found_pyrazole = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pyrazole

        # Check molecule nodes for pyrazole ring
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            if checker.check_ring("pyrazole", mol_smiles):
                found_pyrazole = True
                print(f"Detected pyrazole scaffold in molecule: {mol_smiles}")

        # Check reaction nodes for pyrazole formation
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains pyrazole but reactants don't
                product_has_pyrazole = checker.check_ring("pyrazole", product)

                if product_has_pyrazole:
                    reactants_have_pyrazole = any(
                        checker.check_ring("pyrazole", reactant) for reactant in reactants
                    )

                    if not reactants_have_pyrazole:
                        found_pyrazole = True
                        print(f"Detected pyrazole formation in reaction: {rsmi}")

                # Also check if this is a pyrazole-forming reaction type
                if checker.check_reaction("pyrazole", rsmi) or checker.check_reaction(
                    "{pyrazole}", rsmi
                ):
                    found_pyrazole = True
                    print(f"Detected pyrazole-forming reaction: {rsmi}")
            except Exception as e:
                print(f"Error analyzing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Pyrazole scaffold detected: {found_pyrazole}")
    return found_pyrazole
