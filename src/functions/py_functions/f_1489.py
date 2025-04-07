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
    This function detects a linear synthesis strategy involving a series of
    functional group interconversions without convergent steps.
    """
    # Track the number of functional group transformations
    fg_transformations = 0
    convergent_steps = 0

    # List of functional groups to check for transformations
    fg_types = [
        "Acyl halide",
        "Aldehyde",
        "Ketone",
        "Carboxylic acid",
        "Ester",
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Primary amide",
        "Secondary amide",
        "Tertiary amide",
        "Nitrile",
        "Nitro group",
        "Azide",
        "Alkyne",
        "Alkene",
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Aromatic halide",
        "Phenol",
        "Ether",
        "Thiol",
        "Sulfide",
        "Mesylate",
        "Tosylate",
        "Triflate",
        "Phosphate ester",
    ]

    # Common reagents and solvents to ignore
    common_reagents = [
        "ccoc",
        "cs(=o)",
        "b",
        "na",
        "ccn(cc)cc",
        "clccl",
        "o=p(cl)(cl)cl",
        "o=s(=o)(o)o",
        "cs(=o)(=o)o",
        "cl",
        "o",
    ]

    def is_fg_transformation(reactant_smiles, product_smiles):
        """Check if the reaction involves a functional group transformation"""
        # Check which functional groups are in reactants and product
        reactant_fgs = set()
        product_fgs = set()

        # Clean SMILES strings to handle atom mapping
        clean_product = Chem.MolToSmiles(Chem.MolFromSmiles(product_smiles))

        for fg in fg_types:
            # Check reactants
            for r in reactant_smiles:
                try:
                    clean_r = Chem.MolToSmiles(Chem.MolFromSmiles(r))
                    if checker.check_fg(fg, clean_r):
                        reactant_fgs.add(fg)
                except:
                    continue

            # Check product
            if checker.check_fg(fg, clean_product):
                product_fgs.add(fg)

        # If there's a difference in functional groups, it's a transformation
        return reactant_fgs != product_fgs

    def is_significant_reactant(smiles):
        """Determine if a reactant is significant (not a small reagent)"""
        try:
            # Check if it's a common reagent by SMILES pattern
            if any(reagent in smiles.lower() for reagent in common_reagents):
                return False

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False

            # Count heavy atoms (non-hydrogen)
            heavy_atom_count = mol.GetNumHeavyAtoms()

            # Check for common reagent functional groups
            is_reagent = any(
                checker.check_fg(reagent, smiles)
                for reagent in [
                    "Phosphate ester",
                    "Sulfonamide",
                    "Triflate",
                    "Mesylate",
                    "Tosylate",
                ]
            )

            # Molecules with more than 4 heavy atoms that aren't reagents are considered significant
            return heavy_atom_count > 4 and not is_reagent
        except:
            return False

    def dfs_traverse(node, depth=0):
        nonlocal fg_transformations, convergent_steps

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                reactants = reactants_part.split(".")

                # Count number of significant reactants (excluding small reagents)
                significant_reactants = [
                    r for r in reactants if is_significant_reactant(r)
                ]

                print(f"Reaction at depth {depth}: {rsmi}")
                print(f"Significant reactants: {significant_reactants}")

                # If more than one significant reactant, it's a convergent step
                if len(significant_reactants) > 1:
                    print(f"Convergent step detected at depth {depth}: {rsmi}")
                    convergent_steps += 1

                # Check if this is a functional group transformation
                if is_fg_transformation(reactants, product_part):
                    print(f"FG transformation detected at depth {depth}: {rsmi}")
                    fg_transformations += 1
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if we have multiple FG transformations and no convergent steps
    strategy_detected = fg_transformations >= 3 and convergent_steps == 0
    print(
        f"Linear functional group interconversion strategy detected: {strategy_detected}"
    )
    print(
        f"FG transformations: {fg_transformations}, Convergent steps: {convergent_steps}"
    )
    return strategy_detected
