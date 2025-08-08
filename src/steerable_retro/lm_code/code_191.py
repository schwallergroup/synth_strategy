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
    Detects if a halogenated aromatic system is maintained throughout the synthesis.

    This function checks if at least one halogenated aromatic system is present
    in every molecule node along the main synthetic pathway (excluding reagents).
    """
    # Track if we've found any molecule in the main pathway without a halogenated aromatic
    all_have_halogenated_aromatic = True

    def dfs_traverse(node, depth=0, is_main_path=True):
        nonlocal all_have_halogenated_aromatic

        # Check molecule nodes that are on the main synthetic pathway
        if node["type"] == "mol" and "smiles" in node:
            # Skip checking in-stock starting materials
            if node.get("in_stock", False):
                print(f"Skipping in-stock material at depth {depth}: {node['smiles']}")
                return

            # Only check molecules on the main synthetic pathway
            if is_main_path:
                # Check for aromatic halide using the checker function
                has_aromatic_halide = checker.check_fg("Aromatic halide", node["smiles"])

                if not has_aromatic_halide:
                    print(
                        f"Molecule on main path at depth {depth} does not contain halogenated aromatic: {node['smiles']}"
                    )
                    all_have_halogenated_aromatic = False
                else:
                    print(
                        f"Found aromatic halide in molecule on main path at depth {depth}: {node['smiles']}"
                    )
            else:
                print(f"Skipping reagent/auxiliary molecule at depth {depth}: {node['smiles']}")

        # Process reaction nodes and their children
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Process children
                for child in node.get("children", []):
                    # If this child is a molecule that matches the product, it's on the main path
                    if child["type"] == "mol" and child["smiles"] == product:
                        dfs_traverse(child, depth + 1, True)
                    else:
                        # This is a reagent or auxiliary molecule
                        dfs_traverse(child, depth + 1, False)
            except Exception as e:
                print(f"Error processing reaction node: {e}")
                # If we can't determine the main path, check all children
                for child in node.get("children", []):
                    dfs_traverse(child, depth + 1, is_main_path)
        else:
            # For non-reaction nodes or reaction nodes without proper metadata
            # Process all children with the current path status
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1, is_main_path)

    # Start traversal from the root
    dfs_traverse(route)

    return all_have_halogenated_aromatic
