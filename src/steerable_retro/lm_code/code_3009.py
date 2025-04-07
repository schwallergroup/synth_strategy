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
    This function detects if a synthetic route employs a convergent synthesis strategy
    by joining two complex fragments in a late-stage reaction.
    """
    found_convergent = False

    def calculate_complexity(mol_smiles):
        """Calculate molecular complexity based on atoms, rings, and functional groups"""
        if not mol_smiles or not mol_smiles.strip():
            return 0

        try:
            mol = Chem.MolFromSmiles(mol_smiles)
            if not mol:
                return 0

            # Base complexity: heavy atom count
            complexity = mol.GetNumHeavyAtoms()

            # Add complexity for rings
            ring_info = mol.GetRingInfo()
            complexity += len(ring_info.AtomRings()) * 2

            # Add complexity for functional groups
            fg_list = [
                "Ester",
                "Amide",
                "Alcohol",
                "Amine",
                "Carboxylic acid",
                "Nitrile",
                "Aromatic halide",
                "Aldehyde",
                "Ketone",
                "Alkyne",
                "Alkene",
            ]

            for fg in fg_list:
                if checker.check_fg(fg, mol_smiles):
                    complexity += 2

            return complexity
        except Exception as e:
            print(f"Error calculating complexity: {e}")
            return 0

    def is_coupling_or_joining_reaction(rxn_smiles):
        """Check if the reaction is a coupling or joining reaction commonly used in convergent synthesis"""
        coupling_reactions = [
            "Suzuki coupling",
            "Negishi coupling",
            "Stille reaction",
            "Heck reaction",
            "Sonogashira",
            "Buchwald-Hartwig",
            "N-arylation",
            "Ullmann condensation",
            "Wittig reaction",
            "Diels-Alder",
            "Michael addition",
            "Aldol condensation",
            "Mitsunobu reaction",
            "Huisgen",
            "Ugi reaction",
            "Schotten-Baumann",
        ]

        for rxn_type in coupling_reactions:
            if checker.check_reaction(rxn_type, rxn_smiles):
                print(f"Detected coupling/joining reaction: {rxn_type}")
                return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal found_convergent

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Only consider reactions with at least 2 reactants
            if len(reactants) >= 2:
                # Calculate complexity for each reactant
                reactant_complexities = []

                for reactant in reactants:
                    if reactant.strip():
                        complexity = calculate_complexity(reactant)
                        reactant_complexities.append(complexity)

                # Calculate product complexity
                product_complexity = calculate_complexity(product)

                # Check if we have at least 2 substantial fragments being joined
                # Define "substantial" as having complexity score >= 10
                substantial_fragments = [
                    (i, complexity)
                    for i, complexity in enumerate(reactant_complexities)
                    if complexity >= 10
                ]

                # Check for convergent synthesis conditions:
                # 1. At least 2 substantial fragments
                # 2. Late-stage reaction (depth <= 3)
                # 3. Either a known coupling reaction OR product complexity significantly higher than individual reactants
                if (
                    len(substantial_fragments) >= 2
                    and depth <= 3
                    and (  # Late-stage = depth 0, 1, 2, or 3
                        is_coupling_or_joining_reaction(rsmi)
                        or (product_complexity > 0.8 * sum(reactant_complexities))
                    )
                ):

                    print(f"Found convergent synthesis at depth {depth}: {rsmi}")
                    print(f"Reactant complexities: {reactant_complexities}")
                    print(f"Product complexity: {product_complexity}")
                    print(f"Substantial fragments: {substantial_fragments}")

                    found_convergent = True

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Convergent synthesis strategy detected: {found_convergent}")
    return found_convergent
