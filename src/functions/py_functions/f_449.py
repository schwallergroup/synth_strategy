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
    This function detects if the synthetic route preserves a lactam core
    throughout the synthesis.
    """
    # Track if we've found a lactam in the target molecule
    target_has_lactam = False
    # Track reactions where lactam is preserved
    preserved_reactions = 0
    total_reactions_with_lactam = 0

    # Helper function to check if a molecule contains a lactam
    def contains_lactam(mol_smiles):
        try:
            # Check for various lactam structures (cyclic amides)
            # Lactams can be in different ring sizes
            rings = [
                "cyclopropane",
                "cyclobutane",
                "cyclopentane",
                "cyclohexane",
                "cycloheptane",
                "cyclooctane",
            ]

            # Check for secondary amides in rings (lactams)
            if checker.check_fg("Secondary amide", mol_smiles):
                for ring in rings:
                    if checker.check_ring(ring, mol_smiles):
                        # Direct check if the molecule contains a lactam
                        # This is a simplification that avoids the need to process atom indices
                        mol = Chem.MolFromSmiles(mol_smiles)
                        if mol:
                            # If both secondary amide and ring are present, it's likely a lactam
                            return True

            # Also check for tertiary amides in rings (N-substituted lactams)
            if checker.check_fg("Tertiary amide", mol_smiles):
                for ring in rings:
                    if checker.check_ring(ring, mol_smiles):
                        # If both tertiary amide and ring are present, it's likely a lactam
                        mol = Chem.MolFromSmiles(mol_smiles)
                        if mol:
                            return True

            return False
        except Exception as e:
            print(f"Error in contains_lactam: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal target_has_lactam, preserved_reactions, total_reactions_with_lactam

        try:
            # Check if the current molecule contains a lactam
            if node["type"] == "mol":
                mol_smiles = node["smiles"]

                # If this is the target molecule (depth 0), check for lactam
                if depth == 0:
                    if contains_lactam(mol_smiles):
                        print(f"Target molecule contains lactam: {mol_smiles}")
                        target_has_lactam = True
                    else:
                        print(f"Target molecule does not contain lactam: {mol_smiles}")

            # Check reactions for lactam preservation
            elif node["type"] == "reaction":
                if "metadata" not in node or "rsmi" not in node["metadata"]:
                    return

                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    return

                # Extract product and reactants
                parts = rsmi.split(">")
                if len(parts) < 3:  # Ensure valid reaction SMILES
                    print(f"Invalid reaction SMILES format: {rsmi}")
                    return

                reactants_smiles = parts[0].split(".")
                product_smiles = parts[-1]

                # Check if product contains lactam
                product_has_lactam = contains_lactam(product_smiles)

                # Check if any reactant contains lactam
                reactant_has_lactam = any(
                    contains_lactam(reactant) for reactant in reactants_smiles
                )

                # If either reactant or product has lactam, count this as a reaction involving lactam
                if product_has_lactam or reactant_has_lactam:
                    total_reactions_with_lactam += 1

                    # If both reactant and product have lactam, the core is preserved in this reaction
                    if product_has_lactam and reactant_has_lactam:
                        preserved_reactions += 1
                        print(f"Lactam preserved in reaction: {rsmi}")
                    elif product_has_lactam:
                        print(f"Lactam formed in reaction: {rsmi}")
                    elif reactant_has_lactam:
                        print(f"Lactam lost in reaction: {rsmi}")

            # Traverse children
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)

        except Exception as e:
            print(f"Error in dfs_traverse: {e}")

    # Start traversal
    dfs_traverse(route)

    # Lactam is preserved if:
    # 1. Target molecule has a lactam
    # 2. All reactions involving lactams preserve them
    if target_has_lactam and (
        preserved_reactions == total_reactions_with_lactam
        or total_reactions_with_lactam == 0
    ):
        print("Lactam core preserved throughout synthesis")
        return True

    print(
        f"Lactam preservation: Target has lactam: {target_has_lactam}, Preserved reactions: {preserved_reactions}/{total_reactions_with_lactam}"
    )
    return False
