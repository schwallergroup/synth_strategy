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
    This function detects a synthetic strategy involving multiple amide coupling reactions.

    It traverses the synthesis route and identifies various types of amide coupling reactions,
    including direct coupling of carboxylic acids with amines, acylation with acid halides,
    and other amide formation reactions.

    Returns True if at least 2 amide coupling reactions are found.
    """
    # Count amide couplings
    amide_coupling_count = 0

    # List of reaction types that represent amide couplings
    amide_coupling_reactions = [
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Acyl chloride with ammonia to amide",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with secondary amine to amide",
        "Carboxylic acid with primary amine to amide",
        "Ester with ammonia to amide",
        "Ester with primary amine to amide",
        "Ester with secondary amine to amide",
        "Schotten-Baumann_amide",
        "Carboxylic acid to amide conversion",
    ]

    def count_amide_bonds(mol_smiles):
        """Count the number of amide bonds in a molecule"""
        amide_count = 0
        if checker.check_fg("Primary amide", mol_smiles):
            amide_count += len(checker.get_fg_atom_indices("Primary amide", mol_smiles))
        if checker.check_fg("Secondary amide", mol_smiles):
            amide_count += len(checker.get_fg_atom_indices("Secondary amide", mol_smiles))
        if checker.check_fg("Tertiary amide", mol_smiles):
            amide_count += len(checker.get_fg_atom_indices("Tertiary amide", mol_smiles))
        return amide_count

    def count_amide_formations(reactants, product):
        """Count how many new amide bonds are formed in the reaction"""
        # Count amides in reactants
        reactant_amide_count = sum(count_amide_bonds(r) for r in reactants)

        # Count amides in product
        product_amide_count = count_amide_bonds(product)

        # Return the number of new amide bonds formed
        new_amides = max(0, product_amide_count - reactant_amide_count)
        print(
            f"New amide bonds formed: {new_amides} (Product: {product_amide_count}, Reactants: {reactant_amide_count})"
        )
        return new_amides

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_count

        if node["type"] == "reaction":
            # Extract reaction SMILES
            try:
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if not rsmi:
                    return

                # Extract reactants and product
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                amide_couplings_in_reaction = 0

                # Method 1: Check for known amide coupling reaction types
                for reaction_type in amide_coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found amide coupling reaction: {reaction_type} at depth {depth}")
                        print(f"Reaction SMILES: {rsmi}")
                        amide_couplings_in_reaction += 1
                        # Don't break - continue checking for other reaction types

                # Method 2: Check for amide formation by examining functional groups
                # Always perform this check regardless of reaction type identification

                # Check if product contains amide group
                if (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                ):

                    # Check if reactants contain carboxylic acid, acyl halide, or ester
                    has_acid_derivative = any(
                        checker.check_fg("Carboxylic acid", r)
                        or checker.check_fg("Acyl halide", r)
                        or checker.check_fg("Ester", r)
                        for r in reactants
                    )

                    # Check if reactants contain amine
                    has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Aniline", r)
                        for r in reactants
                    )

                    if has_acid_derivative and has_amine:
                        # Count the actual number of amide bonds formed
                        new_amides = count_amide_formations(reactants, product)
                        if new_amides > 0:
                            print(
                                f"Found {new_amides} amide coupling(s) by functional group analysis at depth {depth}"
                            )
                            print(f"Reaction SMILES: {rsmi}")
                            amide_couplings_in_reaction = max(
                                amide_couplings_in_reaction, new_amides
                            )

                # Add the number of amide couplings found in this reaction
                amide_coupling_count += amide_couplings_in_reaction

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Total amide coupling reactions found: {amide_coupling_count}")

    # Return True if multiple amide couplings are found
    return amide_coupling_count >= 2
