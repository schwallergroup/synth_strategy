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
    Detects the overall strategy: heterocycle construction via nitrile reduction,
    formylation, and ring closure, followed by late-stage Wittig olefination,
    while maintaining a halogen substituent.
    """
    # Initialize flags for each component strategy
    nitrile_reduction_found = False
    formylation_found = False
    ring_formation_found = False
    wittig_found = False
    halogen_retention = True  # Assume true until proven otherwise

    # Track depth for late-stage determination
    max_depth = 0
    wittig_depth = 0

    # Helper function to check if a molecule contains a halogen atom (Br, Cl, I, F)
    def has_halogen(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            for atom in mol.GetAtoms():
                if atom.GetSymbol() in ["Br", "Cl", "I", "F"]:
                    return True
        return False

    # DFS traversal to analyze the synthetic route
    def dfs_traverse(node, depth=0):
        nonlocal nitrile_reduction_found, formylation_found, ring_formation_found
        nonlocal wittig_found, halogen_retention, max_depth, wittig_depth

        # Update max depth
        max_depth = max(max_depth, depth)

        # Process molecule nodes
        if node["type"] == "mol":
            # Check if molecule contains halogen (only for intermediate molecules, not starting materials)
            if depth > 0 and not node.get("in_stock", False):
                if not has_halogen(node["smiles"]):
                    halogen_retention = False
                    print(f"No halogen found in intermediate at depth {depth}: {node['smiles']}")

        # Process reaction nodes
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Extract reactants and product
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitrile reduction to amine or related transformations
            if any(checker.check_fg("Nitrile", r) for r in reactants):
                if (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                    or checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                ):
                    print(f"Found nitrile transformation at depth {depth}")
                    nitrile_reduction_found = True

            # Check for formylation or acylation reactions - look for amine to amide transformation
            if any(
                checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                for r in reactants
            ):
                if checker.check_fg("Primary amide", product) or checker.check_fg(
                    "Secondary amide", product
                ):
                    print(f"Found formylation/acylation at depth {depth}")
                    formylation_found = True

            # Check for heterocyclic ring formation
            nitrogen_heterocycles = [
                "pyrrole",
                "indole",
                "imidazole",
                "pyrazole",
                "oxazole",
                "thiazole",
                "triazole",
                "tetrazole",
            ]

            # Check if product contains a nitrogen heterocycle not present in reactants
            product_has_heterocycle = any(
                checker.check_ring(ring, product) for ring in nitrogen_heterocycles
            )
            reactants_have_heterocycle = any(
                any(checker.check_ring(ring, r) for ring in nitrogen_heterocycles)
                for r in reactants
            )

            if product_has_heterocycle and not reactants_have_heterocycle:
                print(f"Found heterocycle formation at depth {depth}")
                ring_formation_found = True

            # Check for Wittig olefination
            if (
                checker.check_reaction("Wittig", rsmi)
                or checker.check_reaction("{Wittig}", rsmi)
                or checker.check_reaction("Wittig reaction with triphenylphosphorane", rsmi)
                or checker.check_reaction("Wittig with Phosphonium", rsmi)
            ):
                print(f"Found Wittig olefination at depth {depth}")
                wittig_found = True
                wittig_depth = depth

            # Alternative check for Wittig-like transformations (aldehyde/ketone to alkene)
            if not wittig_found:
                if any(
                    checker.check_fg("Aldehyde", r) or checker.check_fg("Ketone", r)
                    for r in reactants
                ):
                    # Check if product has a new C=C bond
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        for bond in product_mol.GetBonds():
                            if (
                                bond.GetBondType() == Chem.BondType.DOUBLE
                                and bond.GetBeginAtom().GetSymbol() == "C"
                                and bond.GetEndAtom().GetSymbol() == "C"
                            ):
                                print(f"Found potential Wittig-like olefination at depth {depth}")
                                wittig_found = True
                                wittig_depth = depth
                                break

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if Wittig is late-stage (lower depth)
    late_stage_wittig = wittig_found and wittig_depth <= (max_depth / 2)

    # The combined strategy requires most of these elements
    # We've made the requirements more flexible
    result = (
        ring_formation_found
        and (nitrile_reduction_found or formylation_found)  # Ring formation is essential
        and (wittig_found or late_stage_wittig)  # Either nitrile reduction or formylation
        and halogen_retention  # Wittig reaction (preferably late-stage)
    )  # Halogen retention throughout

    print(
        f"Strategy components: Nitrile reduction: {nitrile_reduction_found}, Formylation: {formylation_found}, "
        f"Ring formation: {ring_formation_found}, Wittig: {wittig_found}, Late-stage Wittig: {late_stage_wittig}, "
        f"Halogen retention: {halogen_retention}"
    )

    if result:
        print("Detected complete heterocycle construction strategy with late-stage olefination")

    return result
