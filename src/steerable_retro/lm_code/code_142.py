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
    This function detects late-stage Sonogashira coupling to introduce an aryl group.
    """
    sonogashira_found = False

    def dfs_traverse(node, depth=0):
        nonlocal sonogashira_found

        if node["type"] == "reaction" and depth <= 2:  # Consider depths 0, 1, and 2 as late-stage
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check for Sonogashira coupling using the checker function
                sonogashira_reaction_types = [
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

                for reaction_type in sonogashira_reaction_types:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found late-stage Sonogashira coupling: {reaction_type}")

                        # Verify that an aryl group is introduced via Sonogashira coupling
                        try:
                            reactants = rsmi.split(">")[0].split(".")
                            product = rsmi.split(">")[-1]

                            # Look for both aryl halide/triflate and alkyne components
                            aryl_reactant = None
                            alkyne_reactant = None

                            for reactant in reactants:
                                if checker.check_fg("Aromatic halide", reactant):
                                    print(f"Found aromatic halide reactant: {reactant}")
                                    aryl_reactant = reactant
                                elif checker.check_fg("Triflate", reactant):
                                    print(f"Found triflate reactant: {reactant}")
                                    aryl_reactant = reactant

                                if checker.check_fg("Alkyne", reactant):
                                    print(f"Found alkyne reactant: {reactant}")
                                    alkyne_reactant = reactant

                            # Verify both components are present
                            if aryl_reactant is not None and alkyne_reactant is not None:
                                print("Confirmed aryl group introduction in Sonogashira coupling")
                                sonogashira_found = True
                            else:
                                print("Sonogashira coupling found but missing required components")
                        except Exception as e:
                            print(f"Error analyzing reactants and products: {e}")

                # If no reaction type matched, try to analyze the reaction manually
                if not sonogashira_found:
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check for iodoarene or bromoarene in reactants
                        aryl_halide_present = False
                        alkyne_present = False

                        for reactant in reactants:
                            # Check for aromatic halide
                            if "I" in reactant and checker.check_fg("Aromatic halide", reactant):
                                print(f"Found iodoarene reactant: {reactant}")
                                aryl_halide_present = True
                            elif "Br" in reactant and checker.check_fg("Aromatic halide", reactant):
                                print(f"Found bromoarene reactant: {reactant}")
                                aryl_halide_present = True

                            # Check for terminal alkyne
                            if checker.check_fg("Alkyne", reactant):
                                print(f"Found alkyne reactant: {reactant}")
                                alkyne_present = True

                        # Check if product has C≡C-Ar structure (no longer terminal)
                        if (
                            aryl_halide_present
                            and alkyne_present
                            and checker.check_fg("Alkyne", product)
                        ):
                            print("Manual analysis: Found Sonogashira coupling pattern")
                            sonogashira_found = True
                    except Exception as e:
                        print(f"Error in manual analysis: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return sonogashira_found
