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
    Detects if the synthetic route includes a nitration step (introduction of nitro group).
    """
    nitration_found = False

    def dfs_traverse(node):
        nonlocal nitration_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # Check if this is a known nitration reaction type
                nitration_reaction_types = [
                    "Aromatic nitration with HNO3",
                    "Aromatic nitration with NO3 salt",
                    "Aromatic nitration with NO2 salt",
                    "Aromatic nitration with alkyl NO2",
                    "Non-aromatic nitration with HNO3",
                ]

                for rxn_type in nitration_reaction_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Nitration detected: {rxn_type}")
                        nitration_found = True
                        return

                # If no specific nitration reaction type was found, check for nitro group introduction
                has_nitro_in_product = checker.check_fg("Nitro group", product)

                # Check if any reactant has a nitro group
                has_nitro_in_reactants = any(checker.check_fg("Nitro group", r) for r in reactants)

                if has_nitro_in_product and not has_nitro_in_reactants:
                    print("Nitration detected: Nitro group introduced")
                    nitration_found = True
            except Exception as e:
                print(f"Error processing reaction SMILES for nitration detection: {e}")

        # Continue DFS traversal if nitration not found yet
        if not nitration_found:
            for child in node.get("children", []):
                dfs_traverse(child)

    dfs_traverse(route)
    return nitration_found
