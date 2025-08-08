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
    This function detects if the synthetic route involves isoxazole ring formation.
    """
    isoxazole_formation_detected = False

    def dfs_traverse(node):
        nonlocal isoxazole_formation_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if isoxazole is formed in this reaction
                product_has_isoxazole = checker.check_ring("isoxazole", product_smiles)

                # Check if any reactant contains isoxazole
                reactants_with_isoxazole = any(
                    checker.check_ring("isoxazole", r) for r in reactants_smiles
                )

                # Check for relevant reaction types that form isoxazoles
                is_cycloaddition = (
                    checker.check_reaction("Huisgen 1,3 dipolar cycloaddition", rsmi)
                    or checker.check_reaction(
                        "Huisgen alkyne-azide 1,3 dipolar cycloaddition", rsmi
                    )
                    or checker.check_reaction(
                        "Huisgen alkene-azide 1,3 dipolar cycloaddition", rsmi
                    )
                    or checker.check_reaction("[3+2]-cycloaddition of hydrazone and alkyne", rsmi)
                    or checker.check_reaction("[3+2]-cycloaddition of hydrazone and alkene", rsmi)
                    or checker.check_reaction("[3+2]-cycloaddition of diazoalkane and alkyne", rsmi)
                    or checker.check_reaction("[3+2]-cycloaddition of diazoalkane and alkene", rsmi)
                )

                # Isoxazole is formed if it's in the product but not in any reactant
                # and the reaction is a type that can form isoxazoles
                if product_has_isoxazole and not reactants_with_isoxazole:
                    print(f"Isoxazole ring detected in product but not in reactants: {rsmi}")
                    if is_cycloaddition:
                        print(f"Cycloaddition reaction detected that forms isoxazole: {rsmi}")
                        isoxazole_formation_detected = True
                    else:
                        # Check if it's any reaction that forms isoxazole
                        print(f"Checking if reaction forms isoxazole: {rsmi}")
                        isoxazole_formation_detected = True
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return isoxazole_formation_detected
