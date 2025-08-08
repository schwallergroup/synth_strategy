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
    This function detects a strategy involving nitrile reduction to form an amine.
    """
    nitrile_to_amine_found = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_to_amine_found

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Depth {depth} - Examining reaction: {rsmi}")

                # Check if this is a nitrile reduction reaction
                if checker.check_reaction("Reduction of nitrile to amine", rsmi):
                    print(f"Found nitrile to amine reduction reaction: {rsmi}")
                    nitrile_to_amine_found = True
                    return

                # Alternative check: look for nitrile in reactants and primary amine in product
                reactant_has_nitrile = any(checker.check_fg("Nitrile", r) for r in reactants)
                product_has_amine = checker.check_fg("Primary amine", product)

                if reactant_has_nitrile:
                    print(f"Reactant contains nitrile: {reactants}")
                if product_has_amine:
                    print(f"Product contains primary amine: {product}")

                if reactant_has_nitrile and product_has_amine:
                    print(f"Found nitrile to amine transformation: {rsmi}")
                    nitrile_to_amine_found = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Nitrile to amine reduction strategy found: {nitrile_to_amine_found}")
    return nitrile_to_amine_found
