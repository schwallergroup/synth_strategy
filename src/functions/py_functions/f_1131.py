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
    This function detects a functional group interconversion sequence:
    methyl → aldehyde → hydrazone
    """
    # Track the sequence of transformations
    transformation_sequence = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check for functional groups in the molecule
            if checker.check_fg("Hydrazone", mol_smiles):
                transformation_sequence.append(("hydrazone", depth))
                print(f"Detected hydrazone at depth {depth}: {mol_smiles}")

            if checker.check_fg("Aldehyde", mol_smiles):
                transformation_sequence.append(("aldehyde", depth))
                print(f"Detected aldehyde at depth {depth}: {mol_smiles}")

            # For methyl, we need to check if it's attached to an aromatic carbon
            # Since there's no direct "methyl" functional group in the checker
            if mol_smiles and "C" in mol_smiles:
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    methyl_pattern = Chem.MolFromSmarts("c-C")
                    if mol.HasSubstructMatch(methyl_pattern):
                        transformation_sequence.append(("methyl", depth))
                        print(f"Detected methyl group at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for oxidation of methyl to aldehyde
            if any(checker.check_fg("Aldehyde", product) for _ in [1]):
                for reactant in reactants:
                    if Chem.MolFromSmiles(reactant) and Chem.MolFromSmiles(
                        reactant
                    ).HasSubstructMatch(Chem.MolFromSmarts("c-C")):
                        print(f"Detected methyl to aldehyde transformation: {rsmi}")
                        transformation_sequence.append(("methyl_to_aldehyde", depth))

            # Check for aldehyde to hydrazone transformation
            if checker.check_fg("Hydrazone", product):
                for reactant in reactants:
                    if checker.check_fg("Aldehyde", reactant):
                        print(f"Detected aldehyde to hydrazone transformation: {rsmi}")
                        transformation_sequence.append(("aldehyde_to_hydrazone", depth))

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have the required transformations in the correct sequence
    has_methyl = any(t[0] == "methyl" for t in transformation_sequence)
    has_aldehyde = any(t[0] == "aldehyde" for t in transformation_sequence)
    has_hydrazone = any(t[0] == "hydrazone" for t in transformation_sequence)

    # Check for the specific transformations
    has_methyl_to_aldehyde = any(
        t[0] == "methyl_to_aldehyde" for t in transformation_sequence
    )
    has_aldehyde_to_hydrazone = any(
        t[0] == "aldehyde_to_hydrazone" for t in transformation_sequence
    )

    print(f"Transformation sequence: {transformation_sequence}")
    print(
        f"Has methyl: {has_methyl}, Has aldehyde: {has_aldehyde}, Has hydrazone: {has_hydrazone}"
    )
    print(
        f"Has methyl→aldehyde: {has_methyl_to_aldehyde}, Has aldehyde→hydrazone: {has_aldehyde_to_hydrazone}"
    )

    # Verify the sequence is correct by checking depths
    # In retrosynthetic direction, hydrazone should be at lower depth than aldehyde, which should be lower than methyl
    if has_methyl and has_aldehyde and has_hydrazone:
        methyl_depths = [t[1] for t in transformation_sequence if t[0] == "methyl"]
        aldehyde_depths = [t[1] for t in transformation_sequence if t[0] == "aldehyde"]
        hydrazone_depths = [
            t[1] for t in transformation_sequence if t[0] == "hydrazone"
        ]

        if methyl_depths and aldehyde_depths and hydrazone_depths:
            min_methyl_depth = min(methyl_depths)
            min_aldehyde_depth = min(aldehyde_depths)
            min_hydrazone_depth = min(hydrazone_depths)

            # In retrosynthetic direction, hydrazone (final product) should be at lower depth
            # than aldehyde (intermediate), which should be lower than methyl (starting material)
            if min_hydrazone_depth < min_aldehyde_depth < min_methyl_depth:
                print("Correct sequence detected: methyl → aldehyde → hydrazone")
                return True

    # Alternative check: look for the specific transformations
    if has_methyl_to_aldehyde and has_aldehyde_to_hydrazone:
        print(
            "Detected both transformation steps: methyl→aldehyde and aldehyde→hydrazone"
        )
        return True

    return False
