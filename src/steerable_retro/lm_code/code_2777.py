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
    This function detects purine scaffold construction in early stage synthesis.
    """
    early_stage_construction = False

    def dfs_traverse(node, depth=0):
        nonlocal early_stage_construction

        if node["type"] == "reaction" and depth >= 2:  # Early stage (depth >= 2)
            try:
                # Get reaction SMILES
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Extract product and reactants
                product_smiles = rsmi.split(">")[-1]
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Check if product contains purine scaffold
                has_purine_in_product = checker.check_ring("purine", product_smiles)
                print(f"Product contains purine: {has_purine_in_product}")

                # Check if any reactant contains purine scaffold
                reactants_have_purine = False
                for reactant in reactants_smiles:
                    if checker.check_ring("purine", reactant):
                        reactants_have_purine = True
                        print(f"Reactant contains purine: {reactant}")
                        break

                # If product has purine but reactants don't, it's a purine formation
                if has_purine_in_product and not reactants_have_purine:
                    print(f"Detected purine scaffold construction at depth {depth}")
                    early_stage_construction = True

                    # Additional check for specific ring formation reactions, but not required
                    if any(
                        checker.check_reaction(rxn_type, rsmi)
                        for rxn_type in [
                            "Formation of NOS Heterocycles",
                            "{Pictet-Spengler}",
                            "{benzimidazole_derivatives_carboxylic-acid/ester}",
                            "{benzimidazole_derivatives_aldehyde}",
                            "{imidazole}",
                            "{tetrazole_terminal}",
                            "{tetrazole_connect_regioisomere_1}",
                            "{tetrazole_connect_regioisomere_2}",
                        ]
                    ):
                        print(f"Confirmed specific ring formation reaction at depth {depth}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Early stage purine construction detected: {early_stage_construction}")
    return early_stage_construction
