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
    This function detects Fischer indole synthesis as a key ring-forming step.
    Looks for phenylhydrazine reacting with a ketone to form an indole.
    """
    fischer_indole_detected = False

    def dfs_traverse(node):
        nonlocal fischer_indole_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                try:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Analyzing reaction: {rsmi}")

                    # First check if this is a Fischer indole reaction directly
                    if checker.check_reaction("Fischer indole", rsmi):
                        print(f"Fischer indole synthesis detected via reaction check: {rsmi}")
                        fischer_indole_detected = True
                        return

                    # If direct check fails, check for the components
                    has_aryl_hydrazine = False
                    has_ketone = False
                    has_hydrazone = False
                    has_indole_product = False

                    # Check reactants for arylhydrazine, ketone, and hydrazone
                    for reactant in reactants:
                        if not reactant:
                            continue

                        # Check for arylhydrazine (phenylhydrazine or derivatives)
                        if checker.check_fg("Hydrazine", reactant) and any(
                            checker.check_ring(ring, reactant)
                            for ring in ["benzene", "naphthalene", "anthracene"]
                        ):
                            has_aryl_hydrazine = True
                            print(f"Found arylhydrazine in reactant: {reactant}")

                        # Check for ketone
                        if checker.check_fg("Ketone", reactant):
                            has_ketone = True
                            print(f"Found ketone in reactant: {reactant}")

                        # Check for hydrazone (intermediate in Fischer indole synthesis)
                        if checker.check_fg("Hydrazone", reactant) and any(
                            checker.check_ring(ring, reactant)
                            for ring in ["benzene", "naphthalene", "anthracene"]
                        ):
                            has_hydrazone = True
                            print(f"Found aryl hydrazone in reactant: {reactant}")

                    # Check product for indole
                    if product and checker.check_ring("indole", product):
                        has_indole_product = True
                        print(f"Found indole in product: {product}")

                    # Verify components are present for Fischer indole synthesis
                    # Either: arylhydrazine + ketone → indole
                    # Or: aryl hydrazone → indole
                    if has_indole_product and (
                        (has_aryl_hydrazine and has_ketone) or has_hydrazone
                    ):
                        print(f"Fischer indole synthesis detected via component check: {rsmi}")
                        fischer_indole_detected = True

                except Exception as e:
                    print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return fischer_indole_detected
