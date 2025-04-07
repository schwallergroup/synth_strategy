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
    Detects if halogen atoms (F, Cl, Br, I) are preserved throughout the synthesis,
    indicating a strategy where halogens are maintained as key structural elements.
    """
    # Track if we've found a preserved halogen
    preserved_halogen = False

    # Get the target molecule (final product)
    target_mol = Chem.MolFromSmiles(route["smiles"])
    if not target_mol:
        print("Could not parse target molecule")
        return False

    # Check if target has halogens
    halogen_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]")
    if not target_mol.HasSubstructMatch(halogen_pattern):
        print("Target molecule has no halogens")
        return False

    # Find halogen atoms in the target
    target_matches = target_mol.GetSubstructMatches(halogen_pattern)
    if not target_matches:
        print("No halogen matches found in target")
        return False

    print(f"Found {len(target_matches)} halogen atoms in target molecule")

    # Track halogen preservation through the synthesis
    halogen_preserved_count = 0

    def check_halogen_preservation(node, depth=0):
        nonlocal halogen_preserved_count

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if reactants and products have halogens
            reactants_have_halogen = False
            for reactant in reactants_part.split("."):
                if reactant and any(
                    checker.check_fg(fg, reactant)
                    for fg in [
                        "Primary halide",
                        "Secondary halide",
                        "Tertiary halide",
                        "Aromatic halide",
                    ]
                ):
                    reactants_have_halogen = True
                    break

            product_has_halogen = False
            if product_part and any(
                checker.check_fg(fg, product_part)
                for fg in [
                    "Primary halide",
                    "Secondary halide",
                    "Tertiary halide",
                    "Aromatic halide",
                ]
            ):
                product_has_halogen = True

            # If both reactants and product have halogens, this reaction preserves halogens
            if reactants_have_halogen and product_has_halogen:
                print(f"Depth {depth}: Halogen preserved through reaction")
                halogen_preserved_count += 1

            # Check if this is a halogenation reaction (which would add halogens)
            if any(
                checker.check_reaction(rxn, rsmi)
                for rxn in [
                    "Aromatic fluorination",
                    "Aromatic chlorination",
                    "Aromatic bromination",
                    "Aromatic iodination",
                    "Chlorination",
                    "Fluorination",
                    "Iodination",
                    "Bromination",
                ]
            ):
                print(
                    f"Depth {depth}: Detected halogenation reaction - not preservation"
                )
                return False

        # Continue traversing
        for child in node.get("children", []):
            check_halogen_preservation(child, depth + 1)

    # Start traversal from the root
    check_halogen_preservation(route)

    # We need at least 2 reactions preserving halogens to consider it a strategy
    if halogen_preserved_count >= 2:
        print(
            f"Detected halogen preservation throughout synthesis ({halogen_preserved_count} reactions)"
        )
        return True

    # Check if we have at least 3 consecutive molecules with halogens
    consecutive_molecules_with_halogens = 0
    max_consecutive = 0

    def count_consecutive_halogenated_molecules(node):
        nonlocal consecutive_molecules_with_halogens, max_consecutive

        if node["type"] == "mol" and node.get("smiles"):
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(halogen_pattern):
                consecutive_molecules_with_halogens += 1
                max_consecutive = max(
                    max_consecutive, consecutive_molecules_with_halogens
                )
            else:
                consecutive_molecules_with_halogens = 0

        # Continue traversing
        for child in node.get("children", []):
            count_consecutive_halogenated_molecules(child)

    # Reset counter and start traversal again
    consecutive_molecules_with_halogens = 0
    count_consecutive_halogenated_molecules(route)

    if max_consecutive >= 3:
        print(f"Found {max_consecutive} consecutive molecules with halogens")
        return True

    return False
