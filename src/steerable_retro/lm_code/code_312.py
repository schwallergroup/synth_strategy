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
    This function detects an early-stage urea formation from isocyanate.
    Early stage is defined as occurring at depth >= 2 in the synthesis tree.
    """
    urea_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal urea_formation_detected

        if node["type"] == "reaction" and depth >= 2:  # Early stage
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is a urea formation reaction using the checker function
                is_urea_reaction = (
                    checker.check_reaction("Urea synthesis via isocyanate and primary amine", rsmi)
                    or checker.check_reaction(
                        "Urea synthesis via isocyanate and secondary amine", rsmi
                    )
                    or checker.check_reaction("Urea synthesis via isocyanate and diazo", rsmi)
                    or checker.check_reaction("Urea synthesis via isocyanate and sulfonamide", rsmi)
                    or checker.check_reaction("urea", rsmi)
                )

                # Verify functional groups: isocyanate in reactants and urea in product
                isocyanate_in_reactants = False
                for reactant in reactants:
                    try:
                        if checker.check_fg("Isocyanate", reactant):
                            print(f"Found isocyanate in reactant: {reactant}")
                            isocyanate_in_reactants = True
                            break
                    except Exception as e:
                        print(f"Error checking isocyanate in reactant: {e}")
                        continue

                try:
                    urea_in_product = checker.check_fg("Urea", product)
                    if urea_in_product:
                        print(f"Found urea in product: {product}")
                except Exception as e:
                    print(f"Error checking urea in product: {e}")
                    urea_in_product = False

                # Detect urea formation if both conditions are met
                if is_urea_reaction or (isocyanate_in_reactants and urea_in_product):
                    print(f"Detected urea formation at depth {depth}")
                    urea_formation_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Urea formation detected: {urea_formation_detected}")
    return urea_formation_detected
