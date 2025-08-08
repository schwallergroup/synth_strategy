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
    Detects if the synthesis route begins with a halogenated aromatic compound.
    Looks for C-X bond in the earliest synthetic step.
    """
    # First, find all starting materials (in_stock=True)
    starting_materials = []

    def collect_starting_materials(node, depth=0):
        if node["type"] == "mol" and node.get("in_stock", False):
            starting_materials.append((node["smiles"], depth))

        for child in node.get("children", []):
            collect_starting_materials(child, depth + 1)

    collect_starting_materials(route)

    # If no starting materials found, return False
    if not starting_materials:
        print("No starting materials found in the route")
        return False

    # Find the maximum depth (earliest synthetic step)
    max_depth = max(depth for _, depth in starting_materials)
    print(f"Maximum depth found: {max_depth}")

    # Check if any starting material at the maximum depth is a halogenated aromatic compound
    for smiles, depth in starting_materials:
        if depth == max_depth:
            print(f"Checking starting material at depth {depth}: {smiles}")
            if checker.check_fg("Aromatic halide", smiles):
                print(f"Found halogenated aromatic starting material: {smiles}")
                return True

    print("No halogenated aromatic starting materials found at the earliest synthetic step")
    return False
