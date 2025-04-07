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
    This function detects if the synthesis involves a Suzuki coupling reaction
    (aryl halide + boronic acid/ester → C-C bond between aromatic rings).
    """
    suzuki_found = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_found

        # Print node information for debugging
        indent = "  " * depth
        if node["type"] == "mol":
            print(f"{indent}Molecule: {node['smiles']}")
        else:
            print(f"{indent}Reaction node")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"{indent}Checking reaction SMILES: {rsmi}")

            # Use the checker function to directly check for Suzuki coupling reactions
            if checker.check_reaction("Suzuki coupling with boronic acids", rsmi):
                print(f"{indent}Found Suzuki coupling with boronic acids: {rsmi}")
                suzuki_found = True
            elif checker.check_reaction("Suzuki coupling with boronic esters", rsmi):
                print(f"{indent}Found Suzuki coupling with boronic esters: {rsmi}")
                suzuki_found = True
            elif checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi):
                print(f"{indent}Found Suzuki coupling with boronic acids OTf: {rsmi}")
                suzuki_found = True
            elif checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi):
                print(f"{indent}Found Suzuki coupling with boronic esters OTf: {rsmi}")
                suzuki_found = True
            elif checker.check_reaction("Suzuki coupling with sulfonic esters", rsmi):
                print(f"{indent}Found Suzuki coupling with sulfonic esters: {rsmi}")
                suzuki_found = True
            elif checker.check_reaction("{Suzuki}", rsmi):
                print(f"{indent}Found generic Suzuki coupling: {rsmi}")
                suzuki_found = True

            # Fallback: Manual check for Suzuki coupling patterns
            if not suzuki_found:
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for boronic acid/ester in reactants
                    has_boronic = False
                    has_aryl_halide = False

                    for reactant in reactants:
                        if checker.check_fg("Boronic acid", reactant) or checker.check_fg(
                            "Boronic ester", reactant
                        ):
                            print(f"{indent}Found boronic acid/ester in reactant: {reactant}")
                            has_boronic = True

                        # Check for aromatic halide
                        if checker.check_fg("Aromatic halide", reactant):
                            print(f"{indent}Found aromatic halide in reactant: {reactant}")
                            has_aryl_halide = True

                    # If both components are present, it's likely a Suzuki coupling
                    if has_boronic and has_aryl_halide:
                        print(f"{indent}Detected Suzuki coupling pattern manually: {rsmi}")
                        suzuki_found = True
                except Exception as e:
                    print(f"{indent}Error in manual pattern check: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: Suzuki coupling found = {suzuki_found}")
    return suzuki_found
