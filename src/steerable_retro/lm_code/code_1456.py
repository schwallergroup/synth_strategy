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
    This function detects if the synthesis uses a late-stage olefination to introduce a cyanoethylidene group.
    """
    olefination_with_nitrile_detected = False

    # List of olefination reaction types to check
    olefination_reactions = [
        "Wittig",
        "Wittig reaction with triphenylphosphorane",
        "Wittig with Phosphonium",
        "Julia Olefination",
        "O-alkylation of carboxylic acids with diazo compounds",
        "O-alkylation of amides with diazo compounds",
        "Horner-Wadsworth-Emmons",  # Adding HWE reaction which uses phosphonates
    ]

    def dfs_traverse(node, depth=0):
        nonlocal olefination_with_nitrile_detected

        # Print node information for debugging
        if node["type"] == "mol":
            print(f"Examining molecule at depth {depth}: {node['smiles']}")
        elif node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Check if this is a late-stage reaction (depth ≤ 2)
            if depth <= 2:
                # Check if this is an olefination reaction
                is_olefination = False

                # First check if it's a known olefination reaction
                for reaction_type in olefination_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_olefination = True
                        print(f"Olefination reaction detected: {reaction_type}")
                        break

                # If not a known type, check for phosphonate reagents (HWE-type)
                if not is_olefination and "P(=O)" in rsmi.split(">")[0]:
                    print("Potential phosphonate olefination detected")
                    is_olefination = True

                if is_olefination:
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        print(f"Product: {product}")
                        print(f"Reactants: {reactants}")

                        # Check if product contains a nitrile group
                        if checker.check_fg("Nitrile", product):
                            print("Product contains nitrile group")

                            # Check if the product contains an alkene
                            product_mol = Chem.MolFromSmiles(product)
                            alkene_pattern = Chem.MolFromSmarts("C=C")

                            if product_mol and product_mol.HasSubstructMatch(alkene_pattern):
                                print("Product contains alkene")

                                # Check if at least one reactant contains a nitrile group
                                nitrile_in_reactants = False
                                for reactant in reactants:
                                    if checker.check_fg("Nitrile", reactant):
                                        nitrile_in_reactants = True
                                        print(f"Nitrile group found in reactant: {reactant}")
                                        break

                                if nitrile_in_reactants:
                                    # Check if the cyanoethylidene group is newly formed
                                    # Using a more general pattern without specifying bond type
                                    cyanoethylidene_pattern = Chem.MolFromSmarts("C=CC#N")

                                    # Check if cyanoethylidene exists in product
                                    if product_mol.HasSubstructMatch(cyanoethylidene_pattern):
                                        print("Product contains cyanoethylidene group")

                                        # Check if cyanoethylidene exists in any reactant
                                        cyanoethylidene_in_reactants = False
                                        for reactant in reactants:
                                            reactant_mol = Chem.MolFromSmiles(reactant)
                                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                                cyanoethylidene_pattern
                                            ):
                                                cyanoethylidene_in_reactants = True
                                                print(
                                                    f"Cyanoethylidene found in reactant: {reactant}"
                                                )
                                                break

                                        if not cyanoethylidene_in_reactants:
                                            print(
                                                f"Late-stage olefination with nitrile detected at depth {depth}"
                                            )
                                            olefination_with_nitrile_detected = True
                                    else:
                                        print("Product does not contain cyanoethylidene group")
                                        # Try an alternative pattern for cyanoethylidene
                                        alt_pattern = Chem.MolFromSmarts("C=C-C#N")
                                        if product_mol.HasSubstructMatch(alt_pattern):
                                            print(
                                                "Product contains cyanoethylidene group (alternative pattern)"
                                            )
                                            olefination_with_nitrile_detected = True
                    except Exception as e:
                        print(f"Error processing reaction: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(
        f"Late-stage olefination with nitrile strategy detected: {olefination_with_nitrile_detected}"
    )
    return olefination_with_nitrile_detected
