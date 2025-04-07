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
    This function detects a strategy where an aldehyde is formed and then
    transformed via reductive amination to form a C-N bond.
    """
    aldehyde_formation_depth = None
    reductive_amination_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal aldehyde_formation_depth, reductive_amination_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aldehyde formation
                try:
                    product_has_aldehyde = checker.check_fg("Aldehyde", product)

                    if product_has_aldehyde:
                        # Check if any reactant has aldehyde
                        has_aldehyde_in_reactants = any(
                            checker.check_fg("Aldehyde", reactant) for reactant in reactants
                        )

                        if not has_aldehyde_in_reactants:
                            aldehyde_formation_depth = depth
                            print(f"Aldehyde formation detected at depth {depth}")
                except Exception as e:
                    print(f"Error processing aldehyde detection: {e}")

                # Check for reductive amination (aldehyde or ketone to amine transformation)
                try:
                    # Direct check for reductive amination reaction
                    if (
                        checker.check_reaction("Reductive amination with aldehyde", rsmi)
                        or checker.check_reaction("Reductive amination with ketone", rsmi)
                        or checker.check_reaction("Reductive amination with alcohol", rsmi)
                    ):
                        reductive_amination_depth = depth
                        print(f"Reductive amination detected at depth {depth}")
                    else:
                        # Backup check based on functional groups
                        aldehyde_in_reactants = any(
                            checker.check_fg("Aldehyde", reactant) for reactant in reactants
                        )
                        ketone_in_reactants = any(
                            checker.check_fg("Ketone", reactant) for reactant in reactants
                        )

                        # Check for any type of amine in product
                        amine_in_product = (
                            checker.check_fg("Primary amine", product)
                            or checker.check_fg("Secondary amine", product)
                            or checker.check_fg("Tertiary amine", product)
                        )

                        if (aldehyde_in_reactants or ketone_in_reactants) and amine_in_product:
                            reductive_amination_depth = depth
                            print(f"Reductive amination detected at depth {depth}")
                except Exception as e:
                    print(f"Error processing reductive amination detection: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if both aldehyde formation and reductive amination were detected
    # and in the correct sequence (aldehyde formation before reductive amination in synthesis,
    # which means aldehyde formation has higher depth in retrosynthesis)
    if (
        aldehyde_formation_depth is not None
        and reductive_amination_depth is not None
        and aldehyde_formation_depth > reductive_amination_depth
    ):
        print("Aldehyde to amine transformation strategy detected")
        return True
    return False
