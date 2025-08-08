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
    Detects the presence of multiple heterocycles (pyrazine/pyrimidine and indole) in the final product.
    """
    has_pyrazine_or_pyrimidine = False
    has_indole = False

    def dfs_traverse(node, depth=0):
        nonlocal has_pyrazine_or_pyrimidine, has_indole

        if node["type"] == "mol" and depth == 0:  # Final product
            try:
                # Check for pyrazine or pyrimidine
                if checker.check_ring("pyrazine", node["smiles"]):
                    has_pyrazine_or_pyrimidine = True
                    print(f"Pyrazine detected in final product: {node['smiles']}")
                elif checker.check_ring("pyrimidine", node["smiles"]):
                    has_pyrazine_or_pyrimidine = True
                    print(f"Pyrimidine detected in final product: {node['smiles']}")

                if checker.check_ring("indole", node["smiles"]):
                    has_indole = True
                    print(f"Indole detected in final product: {node['smiles']}")

                # Print status for debugging
                if not has_pyrazine_or_pyrimidine:
                    print(f"No pyrazine or pyrimidine detected in final product")
                if not has_indole:
                    print(f"No indole detected in final product")

            except Exception as e:
                print(f"Error processing molecule: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Print final result for debugging
    print(
        f"Final detection results - Pyrazine/Pyrimidine: {has_pyrazine_or_pyrimidine}, Indole: {has_indole}"
    )

    return has_pyrazine_or_pyrimidine and has_indole
