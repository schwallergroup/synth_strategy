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
    This function detects a fragment coupling strategy via N-alkylation
    in the synthetic route.
    """
    n_alkylation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal n_alkylation_found

        print(f"Traversing node at depth {depth}: {node.get('type', 'unknown')}")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction: {rsmi}")

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if we have multiple reactants (indicating potential fragment coupling)
            if len(reactants) >= 2:
                print(f"Found reaction with {len(reactants)} reactants")

                # First check if this is an N-alkylation reaction
                is_n_alkylation = (
                    checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction("Alkylation of amines", rsmi)
                    or checker.check_reaction("Mitsunobu_imide", rsmi)
                    or checker.check_reaction("Mitsunobu_sulfonamide", rsmi)
                    or checker.check_reaction("Mitsunobu aryl ether", rsmi)
                    or checker.check_reaction("Mitsunobu esterification", rsmi)
                    or checker.check_reaction("Mitsunobu aryl ether (intramolecular)", rsmi)
                    or checker.check_reaction("reductive amination", rsmi)
                    or checker.check_reaction("Reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("Reductive amination with ketone", rsmi)
                    or checker.check_reaction("Reductive amination with alcohol", rsmi)
                )

                print(f"Is N-alkylation reaction: {is_n_alkylation}")

                if is_n_alkylation:
                    # Check for amine-containing fragment
                    amine_fragment = None
                    alkylating_fragment = None

                    for reactant in reactants:
                        if (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                            or checker.check_fg("Tertiary amine", reactant)
                        ):
                            amine_fragment = reactant
                            print(f"Found amine fragment: {reactant}")
                        elif (
                            checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                            or checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Tertiary alcohol", reactant)
                            or checker.check_fg("Aldehyde", reactant)
                            or checker.check_fg("Ketone", reactant)
                        ):
                            alkylating_fragment = reactant
                            print(f"Found alkylating fragment: {reactant}")

                    # Check if we found both fragments
                    if amine_fragment and alkylating_fragment:
                        print("Found both amine and alkylating fragments")

                        # Check if both fragments have significant complexity (fragment coupling)
                        amine_mol = Chem.MolFromSmiles(amine_fragment)
                        alkyl_mol = Chem.MolFromSmiles(alkylating_fragment)

                        if amine_mol and alkyl_mol:
                            # Use heavy atom count as a measure of complexity
                            amine_complexity = amine_mol.GetNumHeavyAtoms()
                            alkyl_complexity = alkyl_mol.GetNumHeavyAtoms()

                            print(f"Amine fragment complexity: {amine_complexity} heavy atoms")
                            print(f"Alkylating fragment complexity: {alkyl_complexity} heavy atoms")

                            # Consider it fragment coupling if at least one fragment has significant complexity
                            if amine_complexity >= 6 or alkyl_complexity >= 6:
                                # Verify that the fragments are connected in the product
                                product_mol = Chem.MolFromSmiles(product)

                                if product_mol:
                                    # Check that the product has fewer primary/secondary amines
                                    primary_amine_pattern = Chem.MolFromSmarts("[N;H2]")
                                    secondary_amine_pattern = Chem.MolFromSmarts("[N;H1]")

                                    primary_count_amine = len(
                                        amine_mol.GetSubstructMatches(primary_amine_pattern)
                                    )
                                    secondary_count_amine = len(
                                        amine_mol.GetSubstructMatches(secondary_amine_pattern)
                                    )

                                    primary_count_product = len(
                                        product_mol.GetSubstructMatches(primary_amine_pattern)
                                    )
                                    secondary_count_product = len(
                                        product_mol.GetSubstructMatches(secondary_amine_pattern)
                                    )

                                    print(
                                        f"Primary amines in amine fragment: {primary_count_amine}"
                                    )
                                    print(
                                        f"Secondary amines in amine fragment: {secondary_count_amine}"
                                    )
                                    print(f"Primary amines in product: {primary_count_product}")
                                    print(f"Secondary amines in product: {secondary_count_product}")

                                    # Check if there's a reduction in primary or secondary amines
                                    if (
                                        primary_count_product < primary_count_amine
                                        or (
                                            primary_count_amine > 0
                                            and secondary_count_product > secondary_count_amine
                                        )
                                        or secondary_count_product < secondary_count_amine
                                    ):
                                        print(
                                            f"Found N-alkylation fragment coupling between {amine_fragment} and {alkylating_fragment}"
                                        )
                                        n_alkylation_found = True

                # If we didn't find a direct N-alkylation reaction, check for implicit N-alkylation
                if not n_alkylation_found:
                    print("Checking for implicit N-alkylation")
                    amine_fragment = None
                    alkylating_fragment = None

                    for reactant in reactants:
                        # Check for nitrogen-containing compounds
                        if (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                            or checker.check_fg("Tertiary amine", reactant)
                            or checker.check_fg("Amidinium", reactant)
                            or checker.check_fg("Primary amide", reactant)
                            or checker.check_fg("Secondary amide", reactant)
                            or checker.check_fg("Tertiary amide", reactant)
                            or "N" in reactant
                        ):
                            amine_fragment = reactant
                            print(f"Found potential amine fragment: {reactant}")
                        # Check for potential alkylating agents
                        elif (
                            checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                            or checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Tertiary alcohol", reactant)
                            or checker.check_fg("Aldehyde", reactant)
                            or checker.check_fg("Ketone", reactant)
                            or "Br" in reactant
                            or "Cl" in reactant
                            or "I" in reactant
                        ):
                            alkylating_fragment = reactant
                            print(f"Found potential alkylating fragment: {reactant}")

                    if amine_fragment and alkylating_fragment:
                        print(
                            "Found both potential amine and alkylating fragments in implicit check"
                        )

                        # Check product for reduced NH count or new N-C bonds
                        amine_mol = Chem.MolFromSmiles(amine_fragment)
                        product_mol = Chem.MolFromSmiles(product)
                        alkyl_mol = Chem.MolFromSmiles(alkylating_fragment)

                        if amine_mol and product_mol and alkyl_mol:
                            # Use heavy atom count as a measure of complexity
                            amine_complexity = amine_mol.GetNumHeavyAtoms()
                            alkyl_complexity = alkyl_mol.GetNumHeavyAtoms()

                            print(f"Amine fragment complexity: {amine_complexity} heavy atoms")
                            print(f"Alkylating fragment complexity: {alkyl_complexity} heavy atoms")

                            # Check for primary amine conversion to secondary or tertiary
                            primary_amine_pattern = Chem.MolFromSmarts("[N;H2]")
                            primary_count_reactant = len(
                                amine_mol.GetSubstructMatches(primary_amine_pattern)
                            )
                            primary_count_product = len(
                                product_mol.GetSubstructMatches(primary_amine_pattern)
                            )

                            # Check for secondary amine conversion to tertiary
                            secondary_amine_pattern = Chem.MolFromSmarts("[N;H1]")
                            secondary_count_reactant = len(
                                amine_mol.GetSubstructMatches(secondary_amine_pattern)
                            )
                            secondary_count_product = len(
                                product_mol.GetSubstructMatches(secondary_amine_pattern)
                            )

                            # Check for N-H bonds in reactant and product
                            nh_pattern = Chem.MolFromSmarts("[N;H]")
                            nh_count_reactant = len(amine_mol.GetSubstructMatches(nh_pattern))
                            nh_count_product = len(product_mol.GetSubstructMatches(nh_pattern))

                            print(f"Primary amines in amine fragment: {primary_count_reactant}")
                            print(f"Secondary amines in amine fragment: {secondary_count_reactant}")
                            print(f"Primary amines in product: {primary_count_product}")
                            print(f"Secondary amines in product: {secondary_count_product}")
                            print(f"N-H bonds in amine fragment: {nh_count_reactant}")
                            print(f"N-H bonds in product: {nh_count_product}")

                            # Check for changes in N-H bonds or amine patterns
                            if (
                                primary_count_product < primary_count_reactant
                                or (
                                    primary_count_reactant > 0
                                    and secondary_count_product > secondary_count_reactant
                                )
                                or secondary_count_product < secondary_count_reactant
                                or nh_count_product < nh_count_reactant
                            ):

                                # Consider it fragment coupling if at least one fragment has significant complexity
                                if amine_complexity >= 6 or alkyl_complexity >= 6:
                                    print(
                                        f"Found implicit N-alkylation fragment coupling between {amine_fragment} and {alkylating_fragment}"
                                    )
                                    n_alkylation_found = True

                            # Additional check: look for specific N-alkylation patterns in the product
                            if not n_alkylation_found:
                                # Check if product contains both fragments
                                try:
                                    # Try to find MCS between product and each fragment
                                    amine_mcs = rdFMCS.FindMCS(
                                        [product_mol, amine_mol],
                                        atomCompare=rdFMCS.AtomCompare.CompareElements,
                                        completeRingsOnly=False,
                                    )
                                    alkyl_mcs = rdFMCS.FindMCS(
                                        [product_mol, alkyl_mol],
                                        atomCompare=rdFMCS.AtomCompare.CompareElements,
                                        completeRingsOnly=False,
                                    )

                                    amine_mcs_size = amine_mcs.numAtoms
                                    alkyl_mcs_size = alkyl_mcs.numAtoms

                                    print(f"MCS size with amine fragment: {amine_mcs_size}")
                                    print(f"MCS size with alkylating fragment: {alkyl_mcs_size}")

                                    # If product contains substantial parts of both fragments
                                    if (
                                        amine_mcs_size >= amine_mol.GetNumHeavyAtoms() * 0.7
                                        and alkyl_mcs_size >= alkyl_mol.GetNumHeavyAtoms() * 0.7
                                        and (amine_complexity >= 6 or alkyl_complexity >= 6)
                                    ):
                                        print(f"Found fragment coupling based on MCS analysis")
                                        n_alkylation_found = True
                                except Exception as e:
                                    print(f"Error in MCS analysis: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return n_alkylation_found
