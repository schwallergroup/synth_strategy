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
    Detects Sonogashira coupling (aryl halide + terminal alkyne) in the synthesis route.
    """
    sonogashira_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal sonogashira_detected

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            print(f"Depth {depth}: Checking reaction: {rsmi}")

            # Check for any Sonogashira coupling reaction using the checker
            sonogashira_types = [
                "Sonogashira acetylene_aryl halide",
                "Sonogashira alkyne_aryl halide",
                "Sonogashira acetylene_aryl OTf",
                "Sonogashira alkyne_aryl OTf",
                "Sonogashira acetylene_alkenyl halide",
                "Sonogashira alkyne_alkenyl halide",
                "Sonogashira acetylene_alkenyl OTf",
                "Sonogashira alkyne_alkenyl OTf",
                "Sonogashira acetylene_acyl halide",
                "Sonogashira alkyne_acyl halide",
            ]

            # First try using the checker
            try:
                for rxn_type in sonogashira_types:
                    result = checker.check_reaction(rxn_type, rsmi)
                    print(f"Checking {rxn_type}: {result}")
                    if result:
                        print(f"Sonogashira coupling detected: {rxn_type}")
                        sonogashira_detected = True
                        break
            except Exception as e:
                print(f"Error checking reaction: {e}")

            # If not detected by checker, try manual check
            if not sonogashira_detected:
                try:
                    # Extract reactants and product
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for terminal alkyne in reactants
                    has_terminal_alkyne = False
                    for reactant in reactants:
                        if checker.check_fg("Alkyne", reactant):
                            has_terminal_alkyne = True
                            print(f"Found terminal alkyne in reactant: {reactant}")
                            break

                    # Check for aryl/vinyl halide in reactants
                    has_aryl_halide = False
                    for reactant in reactants:
                        if checker.check_fg("Aromatic halide", reactant) or checker.check_fg(
                            "Alkenyl halide", reactant
                        ):
                            has_aryl_halide = True
                            print(f"Found aryl/vinyl halide in reactant: {reactant}")
                            break

                    # Check for palladium catalyst in reagents
                    has_pd_catalyst = False
                    reagents = rsmi.split(">")[1].split(".")
                    for reagent in reagents:
                        if "Pd" in reagent:
                            has_pd_catalyst = True
                            print(f"Found Pd catalyst in reagents: {reagent}")
                            break

                    # Check for copper co-catalyst in reagents
                    has_cu_catalyst = False
                    for reagent in reagents:
                        if "Cu" in reagent:
                            has_cu_catalyst = True
                            print(f"Found Cu co-catalyst in reagents: {reagent}")
                            break

                    # If we have the key components of a Sonogashira reaction
                    if (
                        has_terminal_alkyne
                        and has_aryl_halide
                        and (has_pd_catalyst or has_cu_catalyst)
                    ):
                        print("Manual detection: Sonogashira coupling pattern detected")
                        sonogashira_detected = True
                except Exception as e:
                    print(f"Error in manual check: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {sonogashira_detected}")
    return sonogashira_detected
