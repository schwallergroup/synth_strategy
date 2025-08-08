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
    This function detects a synthesis route that involves multiple heterocycles
    in the final product.
    """
    # Initialize counters for heterocycles
    heterocycle_count = 0
    heterocycles_found = set()

    # List of heterocycles to check
    heterocycles_to_check = [
        "oxazole",
        "pyrrole",
        "piperazine",
        "imidazole",
        "thiazole",
        "pyridine",
        "furan",
        "thiophene",
        "morpholine",
        "pyrazole",
        "triazole",
        "tetrazole",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_count, heterocycles_found

        if node["type"] == "mol" and depth == 0:  # Final product
            mol_smiles = node["smiles"]
            print(f"Checking final product: {mol_smiles}")

            # Check for various heterocycles
            for heterocycle in heterocycles_to_check:
                if (
                    checker.check_ring(heterocycle, mol_smiles)
                    and heterocycle not in heterocycles_found
                ):
                    heterocycles_found.add(heterocycle)
                    heterocycle_count += 1
                    print(f"{heterocycle.capitalize()} found in final product")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Print summary
    print(f"Found {heterocycle_count} different heterocycles: {', '.join(heterocycles_found)}")

    # Return True if at least 2 different heterocycles were found
    # This threshold can be adjusted based on requirements
    return heterocycle_count >= 2
