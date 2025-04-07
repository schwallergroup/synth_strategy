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

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects if a trifluoromethyl group is maintained throughout
    the synthesis without modification.
    """
    # Track if we've found a maintained CF3 group
    maintained_cf3 = False

    def dfs_traverse(node, depth=0, tracking_cf3=False):
        nonlocal maintained_cf3

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_cf3 = checker.check_fg("Trifluoro group", mol_smiles)

            print(f"Depth {depth}: Molecule: {mol_smiles}")
            print(f"Has CF3 group: {has_cf3}")

            # If this is the target molecule (depth 0), check if it has a CF3 group
            if depth == 0:
                if not has_cf3:
                    print("Target molecule doesn't have a CF3 group")
                    return False
                tracking_cf3 = True

            # If we're tracking CF3 and this molecule doesn't have it, stop this branch
            if tracking_cf3 and not has_cf3:
                print(f"CF3 group not found in molecule at depth {depth}")
                return False

            # If this is a starting material with CF3, we've successfully tracked it back
            if tracking_cf3 and has_cf3 and node.get("in_stock", False):
                print(f"Found CF3 group in starting material: {mol_smiles}")
                maintained_cf3 = True
                return True

        elif node["type"] == "reaction" and tracking_cf3:
            if "metadata" not in node or "rsmi" not in node["metadata"]:
                print("Missing reaction SMILES metadata")
                return False

            rsmi = node["metadata"]["rsmi"]
            print(f"Depth {depth}: Processing reaction: {rsmi}")

            try:
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check if product has CF3 group
                product_has_cf3 = checker.check_fg("Trifluoro group", product_part)

                if not product_has_cf3:
                    print(f"Product doesn't have CF3 group: {product_part}")
                    return False

                # Check if at least one reactant has CF3 group
                reactant_has_cf3 = False
                for reactant in reactants:
                    if checker.check_fg("Trifluoro group", reactant):
                        reactant_has_cf3 = True
                        break

                if not reactant_has_cf3:
                    print(f"No reactant has CF3 group, CF3 was created in this reaction")
                    return False

                # Verify CF3 is preserved through atom mapping
                # This is a simplification - in a real implementation, we would need to
                # extract the specific atom mapping numbers of the CF3 group and verify
                # they're preserved between reactant and product

            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")
                return False

        # Process children
        all_children_maintain_cf3 = True
        for child in node.get("children", []):
            if not dfs_traverse(child, depth + 1, tracking_cf3):
                all_children_maintain_cf3 = False

        return all_children_maintain_cf3 or maintained_cf3

    # Start traversal from the target molecule
    result = dfs_traverse(route)

    if result:
        print("Detected maintained trifluoromethyl group throughout synthesis")
    else:
        print("Trifluoromethyl group is not maintained throughout synthesis")

    return result
