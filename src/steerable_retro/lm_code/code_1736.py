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
    This function detects a synthetic strategy involving azide chemistry,
    specifically the conversion of an amine to an azide group, azide to amine reduction,
    or azide-based click chemistry reactions.
    """
    azide_strategy_found = False

    def dfs_traverse(node, depth=0):
        nonlocal azide_strategy_found

        # Debug information about current node
        indent = "  " * depth
        if node["type"] == "mol":
            print(f"{indent}Examining molecule node: {node['smiles'][:30]}...")

            # Check if this molecule contains an azide group
            if checker.check_fg("Azide", node["smiles"]):
                print(f"{indent}Found molecule with azide group: {node['smiles'][:30]}...")
                # Don't set flag here, as we need to confirm it's part of a strategy

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                print(f"{indent}Examining reaction: {rsmi[:50]}...")

                # Check for various azide-related reactions

                # 1. Amine to azide conversion
                if checker.check_reaction("Amine to azide", rsmi):
                    print(f"{indent}Found amine to azide conversion reaction")
                    azide_strategy_found = True

                # 2. Formation of azides from halogens
                elif checker.check_reaction("Formation of Azides from halogens", rsmi):
                    print(f"{indent}Found formation of azides from halogens")
                    azide_strategy_found = True

                # 3. Formation of azides from boronic acids
                elif checker.check_reaction("Formation of Azides from boronic acids", rsmi):
                    print(f"{indent}Found formation of azides from boronic acids")
                    azide_strategy_found = True

                # 4. Azide reduction (Staudinger)
                elif checker.check_reaction("Azide to amine reduction (Staudinger)", rsmi):
                    print(f"{indent}Found azide reduction via Staudinger reaction")
                    azide_strategy_found = True

                # 5. Huisgen cycloaddition (click chemistry)
                elif any(
                    checker.check_reaction(rxn, rsmi)
                    for rxn in [
                        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                        "Huisgen 1,3 dipolar cycloaddition",
                        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
                        "Huisgen_Cu-catalyzed_1,4-subst",
                        "Huisgen_Ru-catalyzed_1,5_subst",
                        "Huisgen_disubst-alkyne",
                    ]
                ):
                    print(f"{indent}Found Huisgen cycloaddition (click chemistry)")
                    azide_strategy_found = True

                # Fallback checks if specific reaction checks fail
                elif not azide_strategy_found:
                    # Check for azide in product and amine/halide in reactants
                    if checker.check_fg("Azide", product_str):
                        print(f"{indent}Found product with azide group")

                        # Check reactants for potential precursors
                        reactant_smiles_list = reactants_str.split(".")
                        for r_smiles in reactant_smiles_list:
                            if not r_smiles:
                                continue

                            if (
                                checker.check_fg("Primary amine", r_smiles)
                                or checker.check_fg("Secondary amine", r_smiles)
                                or checker.check_fg("Primary halide", r_smiles)
                                or checker.check_fg("Secondary halide", r_smiles)
                                or checker.check_fg("Aromatic halide", r_smiles)
                            ):
                                print(
                                    f"{indent}Found potential azide formation from: {r_smiles[:30]}..."
                                )
                                azide_strategy_found = True
                                break

                    # Check for azide in reactants and products that could result from azide chemistry
                    elif checker.check_fg("Azide", reactants_str):
                        print(f"{indent}Found reactant with azide group")

                        # Check if this might be a click reaction or azide reduction
                        if (
                            checker.check_fg("Primary amine", product_str)
                            or checker.check_fg("Secondary amine", product_str)
                            or checker.check_fg("Triazole", product_str)
                            or checker.check_fg("Tetrazole", product_str)
                        ):
                            print(
                                f"{indent}Found potential azide consumption to: {product_str[:30]}..."
                            )
                            azide_strategy_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    print("Starting traversal of synthetic route...")
    dfs_traverse(route)

    print(f"Azide-based chemistry strategy found: {azide_strategy_found}")
    return azide_strategy_found
