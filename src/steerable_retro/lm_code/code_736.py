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
    This function detects the esterification of a carboxylic acid.
    """
    esterification_found = False

    def dfs_traverse(node):
        nonlocal esterification_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction: {rsmi}")

            # Check for various esterification reaction types
            is_esterification = checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
            is_transesterification = checker.check_reaction("Transesterification", rsmi)
            is_oxidative_esterification = checker.check_reaction(
                "Oxidative esterification of primary alcohols", rsmi
            )

            if is_esterification or is_transesterification or is_oxidative_esterification:
                print(f"Found potential esterification reaction: {rsmi}")

                # Extract reactants and product
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # For standard esterification
                if is_esterification:
                    has_carboxylic_acid = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    )
                    has_alcohol = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                        for r in reactants
                    )
                    has_ester = checker.check_fg("Ester", product)

                    if has_carboxylic_acid and has_ester:
                        print("Confirmed carboxylic acid to ester conversion")
                        esterification_found = True

                # For transesterification
                elif is_transesterification:
                    has_ester_reactant = any(checker.check_fg("Ester", r) for r in reactants)
                    has_alcohol = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                        for r in reactants
                    )
                    has_ester_product = checker.check_fg("Ester", product)

                    if has_ester_reactant and has_alcohol and has_ester_product:
                        print("Confirmed transesterification")
                        esterification_found = True

                # For oxidative esterification
                elif is_oxidative_esterification:
                    has_alcohol = any(checker.check_fg("Primary alcohol", r) for r in reactants)
                    has_ester = checker.check_fg("Ester", product)

                    if has_alcohol and has_ester:
                        print("Confirmed oxidative esterification")
                        esterification_found = True

            # Check for other reactions that might form esters
            elif checker.check_fg("Ester", rsmi.split(">")[-1]):
                print(f"Found reaction producing ester: {rsmi}")

                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acylation reactions that form esters
                if any(checker.check_fg("Acyl halide", r) for r in reactants) and any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    or checker.check_fg("Aromatic alcohol", r)
                    for r in reactants
                ):
                    print("Confirmed acylation to form ester")
                    esterification_found = True

                # Check for O-alkylation reactions that form esters
                elif any(
                    checker.check_fg("Carboxylic acid", r) for r in reactants
                ) and checker.check_fg("Ester", product):
                    print("Confirmed O-alkylation to form ester")
                    esterification_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Esterification found: {esterification_found}")
    return esterification_found
