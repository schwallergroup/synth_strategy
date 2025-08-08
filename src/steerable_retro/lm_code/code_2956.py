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
    Detects a strategy involving the introduction of an alkyne linker to a pyrimidine core
    via a coupling reaction (typically Sonogashira).
    """
    # Track if we found the required reaction
    found_alkyne_coupling = False

    def dfs_traverse(node):
        nonlocal found_alkyne_coupling

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                print("No reaction SMILES found in metadata")
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                print(f"Invalid reaction SMILES format: {rsmi}")
                return

            reactants = parts[0].split(".")
            product = parts[-1]

            # Check if this is a Sonogashira coupling reaction
            is_sonogashira = any(
                [
                    checker.check_reaction("Sonogashira acetylene_aryl halide", rsmi),
                    checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi),
                    checker.check_reaction("Sonogashira acetylene_aryl OTf", rsmi),
                    checker.check_reaction("Sonogashira alkyne_aryl OTf", rsmi),
                ]
            )

            if is_sonogashira:
                print(f"Found potential Sonogashira coupling: {rsmi}")

                # Check for pyrimidine in reactants
                pyrimidine_reactant = None
                alkyne_reactant = None

                for reactant in reactants:
                    try:
                        # Check if reactant has pyrimidine ring
                        if checker.check_ring("pyrimidine", reactant):
                            print(f"Found pyrimidine in reactant: {reactant}")
                            pyrimidine_reactant = reactant

                            # Check if it has a halide (aromatic halide would be present in Sonogashira)
                            if checker.check_fg("Aromatic halide", reactant):
                                print(f"Pyrimidine has aromatic halide")

                        # Check if reactant has alkyne
                        if checker.check_fg("Alkyne", reactant):
                            print(f"Found alkyne in reactant: {reactant}")
                            alkyne_reactant = reactant
                    except Exception as e:
                        print(f"Error checking reactant {reactant}: {e}")
                        continue

                # Check if product has both pyrimidine and alkyne
                if pyrimidine_reactant and alkyne_reactant:
                    try:
                        if checker.check_ring("pyrimidine", product) and checker.check_fg(
                            "Alkyne", product
                        ):
                            print(f"Product contains both pyrimidine and alkyne: {product}")
                            found_alkyne_coupling = True
                    except Exception as e:
                        print(f"Error checking product {product}: {e}")

            # Even if not identified as Sonogashira, check for the structural transformation
            # This handles cases where the reaction might be classified differently
            if not found_alkyne_coupling:
                try:
                    # Check for halogenated pyrimidine in reactants
                    pyrimidine_halide_reactant = None
                    alkyne_reactant = None

                    for reactant in reactants:
                        if checker.check_ring("pyrimidine", reactant) and checker.check_fg(
                            "Aromatic halide", reactant
                        ):
                            print(f"Found halogenated pyrimidine in reactant: {reactant}")
                            pyrimidine_halide_reactant = reactant

                        if checker.check_fg("Alkyne", reactant):
                            print(f"Found alkyne in reactant: {reactant}")
                            alkyne_reactant = reactant

                    # Check if product has pyrimidine with alkyne
                    if (
                        pyrimidine_halide_reactant
                        and alkyne_reactant
                        and checker.check_ring("pyrimidine", product)
                        and checker.check_fg("Alkyne", product)
                    ):
                        print(f"Found structural transformation of alkyne coupling to pyrimidine")
                        found_alkyne_coupling = True
                except Exception as e:
                    print(f"Error in structural check: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_alkyne_coupling
