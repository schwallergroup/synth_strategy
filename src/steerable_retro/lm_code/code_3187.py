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
    This function detects a synthetic strategy that preserves fluorinated aromatic rings
    throughout the synthesis.
    """
    # Track if we've found fluorinated aromatics in the target molecule
    target_has_fluoro_aromatics = False
    # Track if the fluorinated aromatics are preserved throughout synthesis
    preservation_strategy_used = True

    # First check if the target molecule contains fluorinated aromatic rings
    if route["type"] == "mol":
        target_has_fluoro_aromatics = checker.check_fg("Aromatic halide", route["smiles"])
        print(f"Target molecule has fluorinated aromatics: {target_has_fluoro_aromatics}")

        # If target doesn't have fluorinated aromatics, no need to check preservation
        if not target_has_fluoro_aromatics:
            return False

    def dfs_traverse(node, depth=0):
        nonlocal preservation_strategy_used

        # For molecule nodes
        if node["type"] == "mol":
            # Skip further checks if we already determined strategy isn't used
            if not preservation_strategy_used:
                return

            # Check if this molecule has fluorinated aromatic rings
            has_fluoro_aromatics = checker.check_fg("Aromatic halide", node["smiles"])
            print(f"Molecule at depth {depth} has fluorinated aromatics: {has_fluoro_aromatics}")

            # For in-stock molecules, check if they contribute fluorinated aromatics
            if node.get("in_stock", False) and has_fluoro_aromatics:
                print(f"Found fluorinated aromatic in starting material at depth {depth}")

        # For reaction nodes
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Skip further checks if we already determined strategy isn't used
            if not preservation_strategy_used:
                return

            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if any reactant has fluorinated aromatics
                reactants_with_fluoro = any(
                    checker.check_fg("Aromatic halide", r) for r in reactants_smiles
                )
                # Check if product has fluorinated aromatics
                product_has_fluoro = checker.check_fg("Aromatic halide", product_smiles)

                print(f"Reaction at depth {depth}:")
                print(f"  Reactants with fluorinated aromatics: {reactants_with_fluoro}")
                print(f"  Product has fluorinated aromatics: {product_has_fluoro}")

                # If reactants have fluorinated aromatics but product doesn't,
                # or if product has them but no reactant does, then they're not preserved
                if (reactants_with_fluoro and not product_has_fluoro) or (
                    product_has_fluoro and not reactants_with_fluoro
                ):
                    print(f"Fluorinated aromatics not preserved in reaction at depth {depth}")
                    preservation_strategy_used = False
                    return

                # Check if this is a reaction that might modify aromatic fluorine
                # Avoid reactions that typically modify aromatic halides
                if product_has_fluoro and reactants_with_fluoro:
                    if (
                        checker.check_reaction("Aromatic fluorination", rsmi)
                        or checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                        or checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                        )
                    ):
                        print(f"Reaction at depth {depth} might modify aromatic fluorine")
                        preservation_strategy_used = False
                        return
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if target has fluorinated aromatics and they're preserved throughout
    return target_has_fluoro_aromatics and preservation_strategy_used
