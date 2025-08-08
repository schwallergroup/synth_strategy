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
    This function detects the transformation of carboxylic acid to nitrile.
    """
    nitrile_formation_detected = False

    def dfs_traverse(node):
        nonlocal nitrile_formation_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction: {rsmi}")

            # Check for acid to nitrile transformation
            reactants_have_acid = any(
                checker.check_fg("Carboxylic acid", r) for r in reactants if r
            )
            product_has_nitrile = checker.check_fg("Nitrile", product)

            # Check if this is a Schmidt reaction (acid to nitrile)
            is_schmidt_reaction = checker.check_reaction("Schmidt reaction nitrile", rsmi)

            print(
                f"Reactants have acid: {reactants_have_acid}, Product has nitrile: {product_has_nitrile}"
            )
            print(f"Is Schmidt reaction: {is_schmidt_reaction}")

            if reactants_have_acid and product_has_nitrile:
                if is_schmidt_reaction:
                    nitrile_formation_detected = True
                    print(
                        f"Found carboxylic acid to nitrile transformation (Schmidt reaction) at: {rsmi}"
                    )
                else:
                    # Check for other possible reaction types that convert acid to nitrile
                    print("Checking for other acid to nitrile transformations...")
                    # This could be a dehydration of primary amide or other mechanisms
                    # For now, we'll accept any reaction that converts acid to nitrile
                    nitrile_formation_detected = True
                    print(f"Found carboxylic acid to nitrile transformation at: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return nitrile_formation_detected
