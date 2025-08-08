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
    This function detects if Boc protection is maintained throughout the synthesis
    until the final product.

    In retrosynthetic analysis:
    - Depth 0 is the final product (root of the tree)
    - Higher depths are precursors/starting materials

    We need to check if there's at least one continuous path where Boc protection
    is maintained from its introduction until the final product.
    """
    # Track paths where Boc is continuously present
    continuous_boc_paths = []

    # Track current path during traversal
    current_path = []

    def dfs_traverse(node, depth=0, path_has_boc=False):
        nonlocal current_path

        # Add current node to path
        current_path.append((node, depth, path_has_boc))

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if molecule has Boc group
            has_boc = checker.check_fg("Boc", mol_smiles)

            if has_boc:
                print(f"Boc group detected at depth {depth}, SMILES: {mol_smiles[:30]}...")
                path_has_boc = True

            # If this is a leaf node (starting material) and path has Boc
            if not node.get("children") and path_has_boc:
                # Save this path as it maintains Boc protection
                continuous_boc_paths.append(list(current_path))

        elif node["type"] == "reaction":
            # For reaction nodes, check if it's a Boc protection or deprotection
            if "metadata" in node and "rsmi" in node["metadata"]:
                rxn_smiles = node["metadata"]["rsmi"]

                # Check if this is a Boc deprotection reaction
                is_boc_deprotection = checker.check_reaction("Boc amine deprotection", rxn_smiles)

                if is_boc_deprotection:
                    print(f"Boc deprotection reaction detected at depth {depth}")
                    # If we encounter Boc deprotection, this path doesn't maintain Boc protection
                    path_has_boc = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, path_has_boc)

        # Remove current node from path when backtracking
        current_path.pop()

    # Start traversal
    dfs_traverse(route)

    # Check if the final product (depth 0) has Boc
    final_product_has_boc = False
    if route["type"] == "mol" and checker.check_fg("Boc", route["smiles"]):
        final_product_has_boc = True
        print("Final product has Boc group")

    # Check if there's at least one continuous path with Boc protection
    # that reaches the final product
    for path in continuous_boc_paths:
        # Check if this path starts from the final product
        if path[0][1] == 0:  # depth 0
            print(
                "Found a continuous path with Boc protection from final product to a starting material"
            )
            return True

    # If final product doesn't have Boc, that's expected and OK
    if not final_product_has_boc:
        print("Final product doesn't have Boc group (expected if Boc was removed)")

        # Check if Boc was present in the synthesis and properly removed
        boc_present_in_synthesis = any(path_has_boc for _, _, path_has_boc in current_path)
        if boc_present_in_synthesis:
            print("Boc was present in synthesis and properly removed")
            return True

    print("Boc protection was not maintained throughout synthesis")
    return False
