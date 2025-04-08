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
    This function detects a strategy involving multiple heterocycle formations
    (lactam and morpholine) with nitro group reduction in a linear synthesis.
    """
    # Track key features
    has_lactam_formation = False
    has_morpholine_formation = False
    has_epoxide_opening = False
    has_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal has_lactam_formation, has_morpholine_formation, has_epoxide_opening, has_nitro_reduction

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Analyzing reaction: {rsmi}")

            # Check for lactam formation
            if not has_lactam_formation:
                # Check for pyrrolidone (5-membered lactam) or other lactam rings
                lactam_formed = checker.check_ring("pyrrolidone", product_smiles) and not any(
                    checker.check_ring("pyrrolidone", r) for r in reactants_smiles
                )

                # Also check for other lactam formations like 6-membered lactams
                if not lactam_formed:
                    # Check if the reaction forms an amide in a ring
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants_smiles if Chem.MolFromSmiles(r)
                    ]

                    # Look for cyclic amide (lactam) pattern
                    lactam_pattern = Chem.MolFromSmarts("[#6]1[#6][#7][#6](=[#8])[#6]1")
                    if (
                        product_mol
                        and product_mol.HasSubstructMatch(lactam_pattern)
                        and not any(
                            r and r.HasSubstructMatch(lactam_pattern) for r in reactant_mols
                        )
                    ):
                        lactam_formed = True

                if lactam_formed:
                    print("Detected lactam formation")
                    has_lactam_formation = True

            # Check for morpholine formation
            if not has_morpholine_formation:
                # Standard morpholine check
                if checker.check_ring("morpholine", product_smiles) and not any(
                    checker.check_ring("morpholine", r) for r in reactants_smiles
                ):
                    print("Detected morpholine formation")
                    has_morpholine_formation = True
                # Check for morpholine-like structure in the final product
                elif "O1CCN" in product_smiles and not any("O1CCN" in r for r in reactants_smiles):
                    print("Detected morpholine-like structure formation")
                    has_morpholine_formation = True
                # Check if the reaction creates a morpholine-like structure (O-C-C-N in a ring)
                elif not has_morpholine_formation:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants_smiles if Chem.MolFromSmiles(r)
                    ]

                    # Look for morpholine-like pattern (O-C-C-N in a ring)
                    morpholine_pattern = Chem.MolFromSmarts("[#8]1[#6][#6][#7][#6][#6]1")
                    if (
                        product_mol
                        and product_mol.HasSubstructMatch(morpholine_pattern)
                        and not any(
                            r and r.HasSubstructMatch(morpholine_pattern) for r in reactant_mols
                        )
                    ):
                        print("Detected morpholine-like ring formation")
                        has_morpholine_formation = True

            # Check for epoxide opening
            if not has_epoxide_opening:
                if any(
                    checker.check_ring("oxirane", r) for r in reactants_smiles
                ) and not checker.check_ring("oxirane", product_smiles):
                    print("Detected epoxide opening")
                    has_epoxide_opening = True

            # Check for nitro reduction
            if not has_nitro_reduction:
                # Check if any reactant has a nitro group
                if any(checker.check_fg("Nitro group", r) for r in reactants_smiles):
                    # Check if the product has a primary amine that wasn't in the reactants
                    if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                        print("Detected nitro reduction via reaction check")
                        has_nitro_reduction = True
                    elif checker.check_fg("Primary amine", product_smiles) and not any(
                        checker.check_fg("Primary amine", r) for r in reactants_smiles
                    ):
                        print("Detected nitro reduction via functional group change")
                        has_nitro_reduction = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Summary of detected features:")
    print(f"- Lactam formation: {has_lactam_formation}")
    print(f"- Morpholine formation: {has_morpholine_formation}")
    print(f"- Epoxide opening: {has_epoxide_opening}")
    print(f"- Nitro reduction: {has_nitro_reduction}")

    # Count the number of detected features
    features_count = sum(
        [has_lactam_formation, has_morpholine_formation, has_epoxide_opening, has_nitro_reduction]
    )

    # Return True if at least 3 out of 4 features are detected
    return features_count >= 3
