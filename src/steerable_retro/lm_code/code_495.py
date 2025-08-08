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
    Detects if heterocyclic structures (pyridine and pyrazole) are preserved throughout the synthesis.
    """
    # Track molecules at each depth
    molecules_by_depth = {}
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth

        if depth > max_depth:
            max_depth = depth

        if node["type"] == "mol" and "smiles" in node:
            if depth not in molecules_by_depth:
                molecules_by_depth[depth] = []
            molecules_by_depth[depth].append(node["smiles"])

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if at least one heterocycle (pyridine or pyrazole) is present at each depth
    all_depths_preserve_heterocycles = True

    # Find the target molecule (depth 0) and check which heterocycles it contains
    target_has_pyridine = False
    target_has_pyrazole = False

    if 0 in molecules_by_depth:
        for smiles in molecules_by_depth[0]:
            if checker.check_ring("pyridine", smiles):
                target_has_pyridine = True
                print(f"Target molecule contains pyridine")
            if checker.check_ring("pyrazole", smiles):
                target_has_pyrazole = True
                print(f"Target molecule contains pyrazole")

    # If target doesn't have either heterocycle, no need to check preservation
    if not (target_has_pyridine or target_has_pyrazole):
        print("Target molecule doesn't contain pyridine or pyrazole")
        return False

    # Check preservation at each depth
    for depth in range(max_depth + 1):
        if depth not in molecules_by_depth:
            continue

        depth_has_pyridine = False
        depth_has_pyrazole = False

        for smiles in molecules_by_depth[depth]:
            if target_has_pyridine and checker.check_ring("pyridine", smiles):
                depth_has_pyridine = True
            if target_has_pyrazole and checker.check_ring("pyrazole", smiles):
                depth_has_pyrazole = True

        # Check if the heterocycles present in the target are preserved at this depth
        preserved = (not target_has_pyridine or depth_has_pyridine) and (
            not target_has_pyrazole or depth_has_pyrazole
        )

        if not preserved:
            print(f"Depth {depth} does not preserve the heterocycles from the target molecule")
            all_depths_preserve_heterocycles = False
            break

    return all_depths_preserve_heterocycles
