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
    This function detects a synthetic strategy involving sequential functional group
    interconversions from bromide to azide to amine.
    """
    # Initialize flags to track the sequence of functional group interconversions
    br_to_n3_detected = False
    n3_to_nh2_detected = False

    print("Starting analysis for br_to_n3_to_nh2 conversion strategy")

    def dfs_traverse(node, depth=0):
        nonlocal br_to_n3_detected, n3_to_nh2_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Analyzing reaction at depth {depth}: {rsmi}")
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for bromide to azide conversion
                # First check if this is a known azide formation reaction
                if checker.check_reaction("Formation of Azides from halogens", rsmi):
                    print(f"Found 'Formation of Azides from halogens' reaction at depth {depth}")
                    # Verify that at least one reactant contains a bromide
                    has_bromide = False
                    for reactant in reactants:
                        if (
                            checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                            or checker.check_fg("Aromatic halide", reactant)
                        ):
                            # Check if it's specifically a bromide
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                for atom in mol.GetAtoms():
                                    if atom.GetAtomicNum() == 35:  # Bromine atomic number
                                        has_bromide = True
                                        print(f"Found bromide in reactant: {reactant}")
                                        break

                    # Verify that the product contains an azide
                    if has_bromide and checker.check_fg("Azide", product):
                        br_to_n3_detected = True
                        print(f"Detected bromide to azide conversion at depth {depth}")

                # If specific reaction check fails, try a more general approach for bromide to azide
                if not br_to_n3_detected:
                    # Check for bromide in reactants
                    has_bromide = False
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            for atom in mol.GetAtoms():
                                if atom.GetAtomicNum() == 35:  # Bromine atomic number
                                    has_bromide = True
                                    print(f"Found bromide in reactant (general check): {reactant}")
                                    break

                    # Check for azide in product
                    if has_bromide and checker.check_fg("Azide", product):
                        # Additional check to ensure it's a conversion (not just presence)
                        # Look for azide in reactants - if not present, it's likely a conversion
                        azide_in_reactants = any(checker.check_fg("Azide", r) for r in reactants)
                        if not azide_in_reactants:
                            br_to_n3_detected = True
                            print(
                                f"Detected bromide to azide conversion (general) at depth {depth}"
                            )

                # Check for azide to amine conversion
                # Check for azide reduction reactions
                if checker.check_reaction("Azide to amine reduction (Staudinger)", rsmi):
                    print(f"Found 'Azide to amine reduction' reaction at depth {depth}")
                    # Check for azide in reactants
                    has_azide = False
                    for reactant in reactants:
                        if checker.check_fg("Azide", reactant):
                            has_azide = True
                            print(f"Found azide in reactant: {reactant}")
                            break

                    # Check for amine in product
                    if has_azide and checker.check_fg("Primary amine", product):
                        n3_to_nh2_detected = True
                        print(f"Detected azide to amine conversion at depth {depth}")

                # If specific reaction check fails, try a more general approach for azide reduction
                if not n3_to_nh2_detected:
                    # Check for azide in reactants
                    has_azide = False
                    for reactant in reactants:
                        if checker.check_fg("Azide", reactant):
                            has_azide = True
                            print(f"Found azide in reactant (general check): {reactant}")
                            break

                    # Check for amine in product
                    if has_azide and checker.check_fg("Primary amine", product):
                        # Look for evidence of reduction (loss of N2)
                        reactant_mols = [
                            Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
                        ]
                        product_mol = Chem.MolFromSmiles(product)

                        if product_mol:
                            # Count nitrogen atoms in reactants and product
                            n_count_reactants = sum(
                                sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
                                for mol in reactant_mols
                            )
                            n_count_product = sum(
                                1 for atom in product_mol.GetAtoms() if atom.GetAtomicNum() == 7
                            )

                            # If nitrogen count decreased and we have an amine, it's likely a reduction
                            if n_count_reactants > n_count_product:
                                n3_to_nh2_detected = True
                                print(
                                    f"Detected azide to amine conversion (general) at depth {depth}"
                                )

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"br_to_n3_detected: {br_to_n3_detected}, n3_to_nh2_detected: {n3_to_nh2_detected}")
    # Return True if both conversions were detected
    return br_to_n3_detected and n3_to_nh2_detected
