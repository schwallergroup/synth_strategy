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
    Detects isoxazole ring formation from oxime or oxime derivative intermediates.
    """
    isoxazole_formation_found = False
    oxime_intermediate_found = False

    def dfs_traverse(node):
        nonlocal isoxazole_formation_found, oxime_intermediate_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for oxime pattern in reactants
                for reactant in reactants:
                    if checker.check_fg("Oxime", reactant):
                        oxime_intermediate_found = True
                        print(f"Found oxime intermediate: {reactant}")

                # Check if isoxazole is in product but not in reactants
                if checker.check_ring("isoxazole", product):
                    reactant_has_isoxazole = False
                    for reactant in reactants:
                        if checker.check_ring("isoxazole", reactant):
                            reactant_has_isoxazole = True
                            break

                    if not reactant_has_isoxazole:
                        # Check if this is a cycloaddition reaction that could form isoxazole
                        if (
                            checker.check_reaction("Huisgen 1,3 dipolar cycloaddition", rsmi)
                            or checker.check_reaction(
                                "[3+2]-cycloaddition of hydrazone and alkyne", rsmi
                            )
                            or checker.check_reaction(
                                "[3+2]-cycloaddition of hydrazone and alkene", rsmi
                            )
                            or checker.check_reaction(
                                "[3+2]-cycloaddition of diazoalkane and alkyne", rsmi
                            )
                            or checker.check_reaction(
                                "[3+2]-cycloaddition of diazoalkane and alkene", rsmi
                            )
                        ):
                            isoxazole_formation_found = True
                            print(f"Found isoxazole formation in product: {product}")

                        # If no specific reaction type matches, check if we have both oxime and isoxazole
                        # This is a fallback for cases where the reaction type isn't explicitly recognized
                        elif oxime_intermediate_found:
                            isoxazole_formation_found = True
                            print(f"Found isoxazole formation from oxime: {product}")

        elif node["type"] == "mol":
            # Check if molecule contains oxime group
            if checker.check_fg("Oxime", node["smiles"]):
                oxime_intermediate_found = True
                print(f"Found oxime in molecule node: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both conditions are met
    return isoxazole_formation_found and oxime_intermediate_found
