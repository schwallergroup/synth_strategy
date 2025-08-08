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
    This function detects the strategy of maintaining a Boc-protected amine throughout the synthesis.
    """
    boc_present = []

    def dfs_traverse(node, depth=0):
        nonlocal boc_present

        if node["type"] == "mol":
            # Check for Boc group in molecule
            if checker.check_fg("Boc", node["smiles"]):
                boc_present.append(depth)
                print(f"Found Boc group in molecule at depth {depth}")

        elif node["type"] == "reaction":
            # Extract product
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Check for Boc group in product
            if checker.check_fg("Boc", product):
                boc_present.append(depth)
                print(f"Found Boc group in reaction product at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Check if Boc is present throughout synthesis (at the final product and at least at 3 different depths)
    strategy_detected = 0 in boc_present and len(boc_present) >= 3
    print(f"Boc-protected amine strategy detected: {strategy_detected}")
    print(f"Boc found at depths: {sorted(boc_present)}")
    return strategy_detected
