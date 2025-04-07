#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects if the synthetic route involves a trifluoromethyl group
    that remains unchanged throughout the synthesis.
    """
    # Track if we found a trifluoromethyl group that persists through the synthesis
    found_persistent_cf3 = False

    def dfs_traverse(node, depth=0):
        nonlocal found_persistent_cf3

        if node["type"] == "mol":
            smiles = node.get("smiles", "")
            if not smiles:
                return False

            # Check if this molecule has a trifluoromethyl group
            has_cf3 = checker.check_fg("Trifluoro group", smiles)

            # If this is the target molecule (depth 0) and it has CF3, we need to track it
            if depth == 0 and has_cf3:
                # Get the atom indices of all CF3 groups in the target molecule
                cf3_indices = checker.get_fg_atom_indices("Trifluoro group", smiles)
                if cf3_indices:
                    found_persistent_cf3 = True
                    print(f"Found CF3 group in target molecule: {smiles}")
                    return True

            return has_cf3

        elif node["type"] == "reaction":
            # For reaction nodes, check if CF3 is preserved from product to reactants
            try:
                if "metadata" in node and "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]
                    product = rsmi.split(">")[-1]
                    reactants = rsmi.split(">")[0].split(".")

                    # Check if product has CF3
                    product_has_cf3 = checker.check_fg("Trifluoro group", product)

                    # If product doesn't have CF3, no need to check reactants
                    if not product_has_cf3:
                        return False

                    # Check if at least one reactant has CF3
                    reactant_has_cf3 = any(
                        checker.check_fg("Trifluoro group", r) for r in reactants
                    )

                    # If product has CF3 but no reactant does, CF3 was created in this step
                    if product_has_cf3 and not reactant_has_cf3:
                        found_persistent_cf3 = False
                        print(f"CF3 group was created in this reaction: {rsmi}")
                        return False

                    return product_has_cf3
            except Exception as e:
                print(f"Error processing reaction: {e}")
                return False

        # Process children recursively
        all_children_have_cf3 = True
        for child in node.get("children", []):
            child_has_cf3 = dfs_traverse(child, depth + 1)
            all_children_have_cf3 = all_children_have_cf3 and child_has_cf3

        return all_children_have_cf3

    # Start traversal from the root
    result = dfs_traverse(route)

    # We need both a persistent CF3 group and confirmation from traversal
    return found_persistent_cf3 and result
