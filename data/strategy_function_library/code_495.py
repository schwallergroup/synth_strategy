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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

from pathlib import Path
root_data = Path(__file__).parent.parent

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


PRESERVED_HETEROCYCLES = ["pyridine", "pyrazole"]

def main(route):
    """
    Checks for the preservation of specific heterocyclic rings, defined in PRESERVED_HETEROCYCLES, throughout a synthesis. It verifies that if the final product contains a heterocycle from the list, that same heterocycle is present in at least one molecule at every preceding step.
    """
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

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # Depth increases only if current node is NOT a reaction
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Identify which of the specified heterocycles are in the target molecule (depth 0)
    target_heterocycles = set()
    if 0 in molecules_by_depth:
        for smiles in molecules_by_depth[0]:
            for hc in PRESERVED_HETEROCYCLES:
                if checker.check_ring(hc, smiles):
                    target_heterocycles.add(hc)

    # If the target has none of the heterocycles of interest, the strategy is not applicable.
    if not target_heterocycles:
        return False

    # Check that for every depth, all heterocycles found in the target are also present.
    for depth in range(max_depth + 1):
        if depth not in molecules_by_depth:
            continue

        depth_heterocycles = set()
        for smiles in molecules_by_depth[depth]:
            for hc in target_heterocycles: # Optimization: only check for heterocycles we care about
                if checker.check_ring(hc, smiles):
                    depth_heterocycles.add(hc)
        
        # Check if all heterocycles from the target are a subset of the ones at this depth.
        if not target_heterocycles.issubset(depth_heterocycles):
            return False

    return True
