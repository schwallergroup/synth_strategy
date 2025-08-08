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
    Detects if the synthesis route preserves a nitrile functional group throughout
    """
    # Track molecules at each depth to see if they contain nitrile
    molecules_with_nitrile = {}
    molecule_depths = set()

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            smiles = node["smiles"]

            # Add this depth to our set of molecule depths
            molecule_depths.add(depth)

            # Check for nitrile group using the checker function
            if checker.check_fg("Nitrile", smiles):
                molecules_with_nitrile[depth] = molecules_with_nitrile.get(depth, 0) + 1

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Molecules with nitrile at each depth: {molecules_with_nitrile}")
    print(f"All molecule depths: {sorted(molecule_depths)}")

    # Check if at least one molecule at each molecule depth has a nitrile
    nitrile_preserved = True
    for depth in molecule_depths:
        if depth not in molecules_with_nitrile or molecules_with_nitrile[depth] == 0:
            print(f"No nitrile found at depth {depth}")
            nitrile_preserved = False
            break

    return nitrile_preserved
