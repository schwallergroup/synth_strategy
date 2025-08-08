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
    Detects if chloroaromatic pattern is maintained throughout the synthesis.
    This means that once a chloroaromatic group appears, it should be preserved
    in subsequent molecules along the synthesis path.
    """
    # Track molecules with chloroaromatic groups
    molecules_with_chloroaromatic = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Skip starting materials (in_stock)
            if node.get("in_stock", False):
                print(f"Skipping starting material: {mol_smiles}")
                return

            # Check for chloroaromatic pattern using the checker
            has_chloroaromatic = checker.check_fg("Aromatic halide", mol_smiles)

            # Store result with depth information
            molecules_with_chloroaromatic[mol_smiles] = {
                "has_chloroaromatic": has_chloroaromatic,
                "depth": depth,
            }

            print(
                f"Molecule at depth {depth}: {mol_smiles}, has chloroaromatic: {has_chloroaromatic}"
            )

        # Traverse children (going deeper in retrosynthesis)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # If no non-starting materials were found, return False
    if not molecules_with_chloroaromatic:
        print("No non-starting materials found in the route")
        return False

    # Check if the final product (depth 0) has a chloroaromatic group
    final_product = None
    for mol_smiles, info in molecules_with_chloroaromatic.items():
        if info["depth"] == 0:
            final_product = mol_smiles
            break

    if final_product is None:
        print("Could not identify final product")
        return False

    if not molecules_with_chloroaromatic[final_product]["has_chloroaromatic"]:
        print(f"Final product does not have a chloroaromatic group: {final_product}")
        return False

    # Check if chloroaromatic pattern is maintained throughout the synthesis
    # This means if a molecule has a chloroaromatic group, all molecules with lower depth
    # (closer to the final product) should also have it
    for mol_smiles, info in molecules_with_chloroaromatic.items():
        if info["has_chloroaromatic"]:
            # Check all molecules with lower depth
            for other_smiles, other_info in molecules_with_chloroaromatic.items():
                if other_info["depth"] < info["depth"]:
                    if not other_info["has_chloroaromatic"]:
                        print(
                            f"Chloroaromatic pattern lost: {mol_smiles} (depth {info['depth']}) has it, but {other_smiles} (depth {other_info['depth']}) does not"
                        )
                        return False

    print("Chloroaromatic pattern is maintained throughout the synthesis")
    return True
