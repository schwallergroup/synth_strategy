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
    This function detects early introduction of a trifluoroethyl group
    that is carried through the synthesis.
    """
    max_depth = 0
    trifluoroethyl_appearances = {}  # Track depths where trifluoroethyl groups appear

    def dfs_traverse(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "mol" and "smiles" in node:
            # Check for trifluoro group using the checker function
            if checker.check_fg("Trifluoro group", node["smiles"]):
                # Verify it's specifically a trifluoroethyl group (CF3CH2-)
                # by checking for the pattern in the SMILES
                mol_smiles = node["smiles"]
                if "CCC(F)(F)F" in mol_smiles or "CC(F)(F)F" in mol_smiles:
                    print(f"Found trifluoroethyl group at depth {depth} in molecule: {mol_smiles}")
                    trifluoroethyl_appearances[depth] = mol_smiles

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if trifluoroethyl was introduced early and carried through
    if trifluoroethyl_appearances and max_depth > 0:
        # Early introduction means it appears in the deeper half of the synthesis
        # (higher depth values = earlier in synthesis in retrosynthetic analysis)
        early_threshold = max_depth / 2

        # Find the earliest appearance of the trifluoroethyl group
        earliest_appearance = max(trifluoroethyl_appearances.keys())

        # Check if it appears early in the synthesis (deeper than threshold)
        early_introduction = earliest_appearance > early_threshold

        # Check if it persists through to the final product (depth 0)
        persists_to_final = 0 in trifluoroethyl_appearances

        print(f"Earliest trifluoroethyl appearance: depth {earliest_appearance}")
        print(f"Early threshold: {early_threshold}")
        print(f"Early introduction: {early_introduction}")
        print(f"Persists to final product: {persists_to_final}")

        # Return True if it's introduced early and persists to the final product
        return early_introduction and persists_to_final

    print("No trifluoroethyl group found or synthesis has zero depth")
    return False
