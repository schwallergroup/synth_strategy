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

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    Detects if the synthesis uses a convergent approach where a heterocycle
    is introduced to a pre-assembled complex core.
    """
    # Track key features
    found_heterocycle_introduction = False
    found_complex_core = False

    # List of heterocycles to check
    heterocycle_rings = [
        "pyrazole",
        "pyrrole",
        "imidazole",
        "oxazole",
        "thiazole",
        "triazole",
        "tetrazole",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "indole",
        "indazole",
    ]

    # List of heterocycle formation reactions
    heterocycle_reactions = [
        "Paal-Knorr pyrrole synthesis",
        "Fischer indole",
        "benzimidazole_derivatives_aldehyde",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "pyrazole",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_heterocycle_introduction, found_complex_core

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            # Check if this is a late-stage reaction (depth <= 1)
            if depth <= 1:
                try:
                    # Parse reactants and product
                    reactants = reactants_str.split(".")
                    product_mol = Chem.MolFromSmiles(product_str)

                    # Check for heterocycle formation reaction
                    for reaction_type in heterocycle_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(
                                f"Found heterocycle formation reaction: {reaction_type}"
                            )
                            found_heterocycle_introduction = True
                            break

                    # If no specific reaction found, check for heterocycle in reactants
                    if not found_heterocycle_introduction:
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                for ring in heterocycle_rings:
                                    if checker.check_ring(ring, reactant):
                                        print(f"Found heterocycle in reactant: {ring}")
                                        found_heterocycle_introduction = True
                                        break

                    # Check for complex core in reactants (at least 2 rings)
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            ring_info = reactant_mol.GetRingInfo()
                            if ring_info.NumRings() >= 2:
                                print(
                                    f"Found complex core with {ring_info.NumRings()} rings: {reactant}"
                                )
                                found_complex_core = True
                                break

                    # Verify that the heterocycle is actually introduced to the complex core
                    # by checking if the product contains both
                    if (
                        found_heterocycle_introduction
                        and found_complex_core
                        and product_mol
                    ):
                        # Check if product has both heterocycle and multiple rings
                        has_heterocycle_in_product = False
                        for ring in heterocycle_rings:
                            if checker.check_ring(ring, product_str):
                                has_heterocycle_in_product = True
                                break

                        ring_info = product_mol.GetRingInfo()
                        has_complex_core_in_product = ring_info.NumRings() >= 2

                        # If product doesn't have both, reset flags
                        if not (
                            has_heterocycle_in_product and has_complex_core_in_product
                        ):
                            found_heterocycle_introduction = False
                            found_complex_core = False

                except Exception as e:
                    print(f"Error in convergent_heterocycle_strategy: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if both conditions are met
    if found_heterocycle_introduction and found_complex_core:
        print("Found convergent heterocycle strategy")
        return True
    return False
