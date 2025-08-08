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
    Detects if the final step (depth 0) is a Suzuki coupling between an aryl halide
    and a boronic ester/acid to form a biaryl C-C bond.
    """
    final_step_is_suzuki = False
    first_reaction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_suzuki, first_reaction_found

        # If we already found the final step, no need to continue
        if first_reaction_found:
            return

        if node["type"] == "reaction":
            # The first reaction node we encounter is the final step
            first_reaction_found = True

            # Extract reactants and product from the reaction SMILES
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            print(f"Checking final step reaction: {rsmi}")

            # Check if this is a Suzuki coupling reaction using the checker function
            is_suzuki = (
                checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
            )

            if is_suzuki:
                final_step_is_suzuki = True
                print("Detected Suzuki coupling in final step")
                return  # Stop traversal once we find the final step

            # If the checker function didn't identify it as Suzuki, do a manual check
            # This is a fallback in case the reaction isn't properly categorized
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            has_aryl_halide = False
            has_boronic = False

            for reactant in reactants_smiles:
                try:
                    # Check for aromatic halide
                    if checker.check_fg("Aromatic halide", reactant):
                        has_aryl_halide = True
                        print(f"Found aromatic halide in reactant: {reactant}")

                    # Check for boronic acid/ester
                    if checker.check_fg("Boronic acid", reactant) or checker.check_fg(
                        "Boronic ester", reactant
                    ):
                        has_boronic = True
                        print(f"Found boronic acid/ester in reactant: {reactant}")
                except Exception as e:
                    print(f"Error checking reactant {reactant}: {e}")
                    continue

            # If we have both required functional groups, it's likely a Suzuki coupling
            if has_aryl_halide and has_boronic:
                final_step_is_suzuki = True
                print("Detected Suzuki coupling in final step (manual check)")
                return  # Stop traversal

        # Continue traversing children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return final_step_is_suzuki
