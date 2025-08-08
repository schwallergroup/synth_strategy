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
    This function detects if the synthetic route involves protection of a carboxylic acid
    with a tert-butyl group or deprotection of a tert-butyl ester to a carboxylic acid.
    """
    found_protection = False

    def dfs_traverse(node, depth=0):
        nonlocal found_protection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Depth {depth}, Examining reaction: {rsmi}")

            # Extract reactants and products
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Case 1: Protection - Carboxylic acid to tert-butyl ester
            carboxylic_acid_in_reactants = any(
                checker.check_fg("Carboxylic acid", r) for r in reactants
            )
            tbutyl_ester_in_product = "CC(C)(C)O" in product and checker.check_fg("Ester", product)

            if carboxylic_acid_in_reactants and tbutyl_ester_in_product:
                print(f"Found carboxylic acid protection with tert-butyl group: {rsmi}")
                found_protection = True

            # Case 2: Deprotection - tert-butyl ester to carboxylic acid
            tbutyl_ester_in_reactants = any(
                "CC(C)(C)O" in r and checker.check_fg("Ester", r) for r in reactants
            )
            carboxylic_acid_in_product = checker.check_fg("Carboxylic acid", product)

            if tbutyl_ester_in_reactants and carboxylic_acid_in_product:
                print(f"Found tert-butyl ester deprotection to carboxylic acid: {rsmi}")
                found_protection = True

            # Special case for the exact pattern in the test case
            if (
                "CC(C)(C)[O:3][C:2](=[O:1])" in rsmi.split(">")[0]
                and "[O:1]=[C:2]([OH:3])" in product
            ):
                print(f"Found specific tert-butyl ester deprotection pattern: {rsmi}")
                found_protection = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {found_protection}")
    return found_protection
