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
    Detects if the synthesis includes construction of a heterocyclic system.
    """
    heterocycle_formation_detected = False

    # List of heterocyclic rings to check
    heterocycles = [
        "indazole",
        "furan",
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "triazole",
        "tetrazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "thiophene",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
    ]

    def dfs_traverse(node):
        nonlocal heterocycle_formation_detected

        if node["type"] == "reaction":
            try:
                # Get reaction SMILES
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Split reactants
                reactants = reactants_part.split(".")

                # Check for heterocycle formation
                for heterocycle in heterocycles:
                    # Check if heterocycle is in product
                    if checker.check_ring(heterocycle, product):
                        # Check if heterocycle is not in any reactant
                        reactant_has_heterocycle = any(
                            checker.check_ring(heterocycle, r) for r in reactants
                        )

                        if not reactant_has_heterocycle:
                            print(f"Detected heterocycle ({heterocycle}) formation")
                            heterocycle_formation_detected = True
                            break
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return heterocycle_formation_detected
