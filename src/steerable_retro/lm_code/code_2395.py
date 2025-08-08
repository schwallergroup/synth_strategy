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
    This function detects a linear synthesis strategy involving heterocyclic compounds
    (specifically pyrazole and piperidine rings).
    """
    # Track heterocycles
    contains_pyrazole = False
    contains_piperidine = False
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal contains_pyrazole, contains_piperidine, is_linear

        print(f"Checking node at depth {depth}: {node['type']}")

        if node["type"] == "mol":
            # Check for heterocycles in the molecule using the checker functions
            if checker.check_ring("pyrazole", node["smiles"]):
                contains_pyrazole = True
                print(f"Found pyrazole at depth {depth} in molecule: {node['smiles']}")

            if checker.check_ring("piperidine", node["smiles"]):
                contains_piperidine = True
                print(f"Found piperidine at depth {depth} in molecule: {node['smiles']}")

        elif node["type"] == "reaction":
            # Check if reaction has multiple products (non-linear)
            try:
                rsmi = node["metadata"]["rsmi"]
                products = rsmi.split(">")[-1].split(".")
                if len(products) > 1:
                    is_linear = False
                    print(f"Found non-linear reaction at depth {depth}: {rsmi}")
            except KeyError:
                print(f"Warning: No rsmi found in reaction metadata at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting traversal of synthesis route")
    dfs_traverse(route)

    # Print final results for debugging
    print(
        f"Final results: contains_pyrazole={contains_pyrazole}, contains_piperidine={contains_piperidine}, is_linear={is_linear}"
    )

    # Return True if it's a linear synthesis with both heterocycles
    return is_linear and contains_pyrazole and contains_piperidine
