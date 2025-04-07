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
    This function detects if a synthesis route includes a sequence of
    primary amine → azide → carbamate transformations.
    """
    # Track molecules and their transformations with atom mapping and depth
    azide_to_amine_transformations = []  # Store (product_mol, azide_reactant, depth)
    amine_to_carbamate_transformations = (
        []
    )  # Store (amine_reactant, carbamate_product, depth)

    def dfs_traverse(node, depth=0):
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")
            product = rsmi.split(">")[-1]

            print(f"Depth {depth}, Examining reaction: {rsmi}")

            try:
                # Check for azide→amine transformation (azide in reactant, amine in product)
                for reactant in reactants:
                    if checker.check_fg("Azide", reactant):
                        print(f"Found azide in reactant: {reactant}")
                        if checker.check_fg("Primary amine", product):
                            print(f"Found primary amine in product: {product}")
                            # Store the transformation with atom mapping and depth
                            azide_to_amine_transformations.append(
                                (product, reactant, depth)
                            )
                            print(
                                f"Recorded azide→amine transformation at depth {depth}: {reactant} → {product}"
                            )

                # Check for amine→carbamate transformation (amine in reactant, carbamate in product)
                for reactant in reactants:
                    if checker.check_fg("Primary amine", reactant):
                        print(f"Found primary amine in reactant: {reactant}")
                        if checker.check_fg("Carbamic ester", product):
                            print(f"Found carbamate in product: {product}")
                            # Store the transformation with atom mapping and depth
                            amine_to_carbamate_transformations.append(
                                (reactant, product, depth)
                            )
                            print(
                                f"Recorded amine→carbamate transformation at depth {depth}: {reactant} → {product}"
                            )
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children (moving backward in synthesis)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Azide→amine transformations: {len(azide_to_amine_transformations)}")
    print(f"Amine→carbamate transformations: {len(amine_to_carbamate_transformations)}")

    # Check if we have the sequence on the same molecule
    sequence_found = False

    # For each azide→amine transformation
    for amine_product, azide_reactant, azide_depth in azide_to_amine_transformations:
        # For each amine→carbamate transformation
        for (
            amine_reactant,
            carbamate_product,
            carbamate_depth,
        ) in amine_to_carbamate_transformations:
            # Check if the amine product from azide reduction is used as reactant in carbamate formation
            # and ensure the azide→amine happens before amine→carbamate (higher depth means earlier in synthesis)
            if azide_depth > carbamate_depth:
                try:
                    # Extract atom mapping numbers to track the same nitrogen atom
                    amine_product_mol = Chem.MolFromSmiles(amine_product)
                    amine_reactant_mol = Chem.MolFromSmiles(amine_reactant)

                    if amine_product_mol and amine_reactant_mol:
                        # Use MCS to check if they're the same core structure
                        mcs = rdFMCS.FindMCS(
                            [amine_product_mol, amine_reactant_mol],
                            atomCompare=rdFMCS.AtomCompare.CompareElements,
                            bondCompare=rdFMCS.BondCompare.CompareOrder,
                            ringMatchesRingOnly=True,
                            completeRingsOnly=True,
                        )

                        mcs_match_ratio = mcs.numAtoms / min(
                            amine_product_mol.GetNumAtoms(),
                            amine_reactant_mol.GetNumAtoms(),
                        )
                        print(
                            f"MCS match ratio: {mcs_match_ratio:.2f} between {amine_product} and {amine_reactant}"
                        )

                        if mcs_match_ratio >= 0.8:
                            print(
                                f"Found matching molecules in sequence: {amine_product} and {amine_reactant}"
                            )
                            print(
                                f"Sequence found: azide→amine at depth {azide_depth}, amine→carbamate at depth {carbamate_depth}"
                            )
                            sequence_found = True
                            break
                except Exception as e:
                    print(f"Error comparing molecules: {e}")

        if sequence_found:
            break

    # Alternative check: look for a reaction that directly converts azide to carbamate
    for node in route.get("children", []):
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")
            product = rsmi.split(">")[-1]

            for reactant in reactants:
                if checker.check_fg("Azide", reactant) and checker.check_fg(
                    "Carbamic ester", product
                ):
                    print(
                        f"Found direct azide→carbamate transformation: {reactant} → {product}"
                    )
                    sequence_found = True

    print(f"Azide → amine → carbamate sequence detected: {sequence_found}")
    return sequence_found
