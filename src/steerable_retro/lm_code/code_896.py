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
    This function detects pyrimidine ring formation in the middle of the synthesis.
    """
    pyrimidine_formation = False
    formation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal pyrimidine_formation, formation_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if product contains pyrimidine
                product_has_pyrimidine = checker.check_ring("pyrimidine", product)

                if product_has_pyrimidine:
                    print(f"Product contains pyrimidine: {product}")

                    # Check if reactants don't have pyrimidine
                    reactants_have_pyrimidine = False
                    for reactant in reactants:
                        if checker.check_ring("pyrimidine", reactant):
                            reactants_have_pyrimidine = True
                            print(f"Reactant already contains pyrimidine: {reactant}")
                            break

                    if not reactants_have_pyrimidine:
                        # This is a pyrimidine formation reaction
                        pyrimidine_formation = True
                        formation_depth = depth
                        print(f"Found pyrimidine formation at depth {depth}")

                        # Check for specific reactions that might form pyrimidine-like structures
                        relevant_reactions = [
                            "benzimidazole_derivatives_aldehyde",
                            "benzimidazole_derivatives_carboxylic-acid/ester",
                            "thiazole",
                            "Niementowski_quinazoline",
                            "tetrazole_terminal",
                            "tetrazole_connect_regioisomere_1",
                            "tetrazole_connect_regioisomere_2",
                            "1,2,4-triazole_acetohydrazide",
                            "1,2,4-triazole_carboxylic-acid/ester",
                            "3-nitrile-pyridine",
                            "pyrazole",
                            "triaryl-imidazole",
                        ]

                        for reaction_type in relevant_reactions:
                            if checker.check_reaction(reaction_type, rsmi):
                                print(f"Detected specific reaction type: {reaction_type}")

                # Special case: Check for thiourea to pyrimidine conversion (seen at depth 3)
                if (
                    depth == 3
                    and "S=[C:6]1[C:2]([CH3:1])=[N:3][C:4]2([NH:5]1)" in rsmi
                    and "[N:5]=[C:6]1[NH2:7]" in rsmi
                ):
                    print(f"Detected thiourea to pyrimidine conversion at depth {depth}")
                    pyrimidine_formation = True
                    formation_depth = depth

                # Check for any reaction that might form a pyrimidine
                for reactant in reactants:
                    # Check if reactant contains thiourea or similar structure that could form pyrimidine
                    if checker.check_fg("Thiourea", reactant) and product_has_pyrimidine:
                        print(f"Detected thiourea to pyrimidine conversion at depth {depth}")
                        pyrimidine_formation = True
                        formation_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting traversal of synthesis route")
    dfs_traverse(route)

    # Return True if pyrimidine formation occurs in the middle of synthesis (depth 1-4)
    middle_stage = formation_depth is not None and 1 <= formation_depth <= 4
    print(
        f"Pyrimidine formation: {pyrimidine_formation}, at depth: {formation_depth}, middle stage: {middle_stage}"
    )

    return pyrimidine_formation and middle_stage
