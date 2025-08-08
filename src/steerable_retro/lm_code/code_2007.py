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
    This function detects the preservation of a cyano group throughout the synthesis,
    indicating its importance as an electron-withdrawing or directing group.
    """
    # Track cyano group presence at each depth
    cyano_at_depth = {}

    def dfs_traverse(node, current_depth=0):
        # For molecule nodes, check for cyano group
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_cyano = checker.check_fg("Nitrile", mol_smiles)

            if has_cyano:
                cyano_at_depth[current_depth] = True
                print(f"Found cyano group in molecule at depth {current_depth}: {mol_smiles}")

        # For reaction nodes, check if cyano is preserved from reactants to product
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                reactants = reactants_part.split(".")

                # Check for cyano in reactants and product
                reactant_has_cyano = any(checker.check_fg("Nitrile", r) for r in reactants)
                product_has_cyano = checker.check_fg("Nitrile", product)

                # Get depth information
                if "Depth:" in node.get("metadata", {}).get("ID", ""):
                    depth_str = (
                        node.get("metadata", {}).get("ID", "").split("Depth:")[1].strip().split()[0]
                    )
                    try:
                        depth = int(depth_str)
                        current_depth = depth  # Update current depth
                    except ValueError:
                        print(f"Could not parse depth from: {depth_str}")

                # Record cyano presence at this depth
                if product_has_cyano:
                    cyano_at_depth[current_depth] = True
                    print(f"Found cyano group in product at depth {current_depth}: {product}")

                # Check if cyano is preserved in this reaction
                if reactant_has_cyano and product_has_cyano:
                    print(f"Cyano group preserved in reaction at depth {current_depth}")
                elif reactant_has_cyano and not product_has_cyano:
                    print(f"Cyano group lost in reaction at depth {current_depth}")
                elif not reactant_has_cyano and product_has_cyano:
                    print(f"Cyano group introduced in reaction at depth {current_depth}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if cyano group is preserved throughout synthesis
    # We need at least 3 depths with cyano group to consider it preserved
    consecutive_depths = 0
    max_consecutive = 0

    for depth in sorted(cyano_at_depth.keys()):
        if depth - 1 in cyano_at_depth:
            consecutive_depths += 1
        else:
            consecutive_depths = 1

        max_consecutive = max(max_consecutive, consecutive_depths)

    print(f"Maximum consecutive depths with cyano group: {max_consecutive}")
    print(f"Total depths with cyano group: {len(cyano_at_depth)}")

    # Return True if cyano group is found at 3 or more depths
    if len(cyano_at_depth) >= 3:
        print("Cyano group preserved throughout synthesis")
        return True
    return False
