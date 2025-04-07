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
    This function detects if an ester group is preserved throughout the entire synthesis.
    It checks if at least one ester group is maintained without modification through all reaction steps.
    """
    # Track if we've found at least one ester that's preserved throughout
    ester_preserved = False

    def check_ester_preservation(node, depth=0):
        nonlocal ester_preserved

        if node["type"] == "mol":
            # Check if this molecule has an ester
            current_has_ester = checker.check_fg("Ester", node["smiles"])
            print(f"Depth {depth}, Molecule: {node['smiles']}, Has ester: {current_has_ester}")

            # If this is a starting material and has an ester, we've found a preserved path
            if node.get("in_stock", False) and current_has_ester:
                print(f"Found starting material with ester: {node['smiles']}")
                return True

            # If no children or no ester, this path doesn't preserve an ester
            if not current_has_ester or not node.get("children"):
                return False

            # Check if any child path preserves the ester
            for child in node.get("children", []):
                if check_ester_preservation(child, depth + 1):
                    return True

            # If we get here, no child path preserved the ester
            return False

        elif node["type"] == "reaction":
            # For reaction nodes, check if ester is preserved through the reaction
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["rsmi"]
                    # In retrosynthesis, the product in rsmi is our current molecule
                    # and the reactants are what we're breaking down into
                    product = rsmi.split(">")[-1]
                    reactants = rsmi.split(">")[0].split(".")

                    # Check if ester is in both product and at least one reactant
                    product_has_ester = checker.check_fg("Ester", product)
                    reactants_with_ester = [r for r in reactants if checker.check_fg("Ester", r)]

                    print(f"Depth {depth}, Reaction: {rsmi}")
                    print(f"  Product has ester: {product_has_ester}")
                    print(f"  Reactants with ester: {len(reactants_with_ester)}")

                    # If product doesn't have ester or no reactant has ester, this path doesn't preserve
                    if not product_has_ester or not reactants_with_ester:
                        return False

                    # Check if any child path preserves the ester
                    for child in node.get("children", []):
                        if check_ester_preservation(child, depth + 1):
                            return True

                    # If we get here, no child path preserved the ester
                    return False

            except Exception as e:
                print(f"Error processing reaction: {e}")
                return False

            # Process children if we haven't returned yet
            for child in node.get("children", []):
                if check_ester_preservation(child, depth + 1):
                    return True

            return False

    # Start traversal from the root (final product)
    if route["type"] == "mol":
        # Check if final product has an ester
        if checker.check_fg("Ester", route["smiles"]):
            print(f"Final product has ester: {route['smiles']}")
            # Check if this ester can be traced back to a starting material
            ester_preserved = check_ester_preservation(route)
        else:
            print(f"Final product has no ester: {route['smiles']}")
    else:
        # If root is a reaction node, process it normally
        ester_preserved = check_ester_preservation(route)

    return ester_preserved
