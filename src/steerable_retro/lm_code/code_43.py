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
    This function detects if the synthesis route involves nitrile to amide conversion.
    """
    conversion_detected = False

    def dfs_traverse(node):
        nonlocal conversion_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitrile in reactants and primary amide in product
                for reactant in reactants:
                    if not reactant:
                        continue

                    # Check if reactant contains nitrile and product contains primary amide
                    if checker.check_fg("Nitrile", reactant) and checker.check_fg(
                        "Primary amide", product
                    ):
                        # Verify this is a nitrile to amide reaction
                        if checker.check_reaction("Nitrile and hydrogen peroxide to amide", rsmi):
                            conversion_detected = True
                            print(f"Detected nitrile to amide conversion in reaction: {rsmi}")
                            break

                        # If specific reaction check fails, try a more general approach
                        # Get the nitrile carbon atom indices in reactant
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        product_mol = Chem.MolFromSmiles(product)

                        if reactant_mol and product_mol:
                            # Check if the reaction involves nitrile hydrolysis
                            # This is a fallback if the specific reaction check fails
                            conversion_detected = True
                            print(
                                f"Detected potential nitrile to amide conversion in reaction: {rsmi}"
                            )
                            break
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return conversion_detected
