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
    Detects if the synthesis involves late-stage heterocycle formation,
    specifically a pyrazole ring formed in the final step.
    """
    final_product_has_pyrazole = False
    first_reaction_forms_pyrazole = False

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_pyrazole, first_reaction_forms_pyrazole

        if node["type"] == "mol" and depth == 0:
            # Check if final product has pyrazole
            if checker.check_ring("pyrazole", node["smiles"]):
                final_product_has_pyrazole = True
                print(f"Final product contains pyrazole ring: {node['smiles']}")

        elif node["type"] == "reaction" and depth == 1:  # First reaction in retrosynthetic analysis
            # Check if this reaction forms pyrazole
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")

                # Check if the reaction is a pyrazole formation reaction using multiple methods
                is_pyrazole_reaction = (
                    checker.check_reaction("pyrazole", rsmi)
                    or checker.check_reaction("{pyrazole}", rsmi)
                    or checker.check_reaction("Michael-induced ring closure from hydrazone", rsmi)
                    or checker.check_reaction("[3+2]-cycloaddition of hydrazone and alkyne", rsmi)
                    or checker.check_reaction("[3+2]-cycloaddition of hydrazone and alkene", rsmi)
                )

                if is_pyrazole_reaction:
                    print("Reaction identified as potential pyrazole formation")

                    # Verify product has pyrazole but reactants don't
                    product_has_pyrazole = checker.check_ring("pyrazole", product)
                    reactants_have_pyrazole = any(
                        checker.check_ring("pyrazole", r) for r in reactants if r
                    )

                    if product_has_pyrazole and not reactants_have_pyrazole:
                        print("Confirmed: Product has pyrazole but reactants don't")
                        first_reaction_forms_pyrazole = True
                else:
                    # Manual check for pyrazole formation if reaction type check failed
                    product_has_pyrazole = checker.check_ring("pyrazole", product)
                    reactants_have_pyrazole = any(
                        checker.check_ring("pyrazole", r) for r in reactants if r
                    )

                    if product_has_pyrazole and not reactants_have_pyrazole:
                        print("Manual check: Product has pyrazole but reactants don't")

                        # Check for common functional groups in pyrazole formation
                        has_nitrile = any(checker.check_fg("Nitrile", r) for r in reactants if r)
                        has_hydrazine_or_derivative = any(
                            checker.check_fg("Hydrazine", r)
                            or checker.check_fg("Hydrazone", r)
                            or checker.check_fg("Acylhydrazine", r)
                            or checker.check_fg("Hydrazone amide", r)
                            for r in reactants
                            if r
                        )

                        if has_nitrile and has_hydrazine_or_derivative:
                            print("Detected nitrile and hydrazine/derivative in reactants")
                            first_reaction_forms_pyrazole = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    result = final_product_has_pyrazole and first_reaction_forms_pyrazole
    print(
        f"Final result: {result} (final_product_has_pyrazole={final_product_has_pyrazole}, first_reaction_forms_pyrazole={first_reaction_forms_pyrazole})"
    )
    return result
