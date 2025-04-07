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
    This function detects if the synthetic route uses a linear assembly of multiple fragments.
    Linear assembly involves sequential coupling of fragments, where most steps involve
    combining two significant molecular fragments.
    """
    fragment_counts = []
    coupling_reactions = []

    def is_reagent_or_catalyst(smiles):
        """Check if a SMILES string represents a common reagent or catalyst"""
        # Common reagents and catalysts patterns
        reagent_patterns = [
            r"^\[.*\]$",  # Metal catalysts like [Pd]
            r"^[A-Z][a-z]?[+-]?$",  # Simple ions like Na+, Cl-
            r"^O$",
            r"^CO$",
            r"^CCO$",
            r"^CCOC\(C\)=O$",  # Common solvents
            r"^CN\(C\)C=O$",
            r"^O=C\(O\)O$",  # DMF, carboxylic acids
            r"^ClCCl$",
            r"^\[H\]\[H\]$",  # DCM, H2
            r"^CCN\(C\(C\)C\)C\(C\)C$",  # DIPEA
            r"^B$",  # Boron
        ]

        # Check if the SMILES matches any reagent pattern
        for pattern in reagent_patterns:
            if re.match(pattern, smiles):
                return True

        # Check for small molecules (likely reagents)
        mol = Chem.MolFromSmiles(smiles)
        if mol and mol.GetNumAtoms() < 5:
            return True

        return False

    def is_coupling_reaction(rsmi):
        """Check if a reaction is a coupling reaction"""
        # Common coupling reaction types
        coupling_types = [
            "Suzuki",
            "Buchwald-Hartwig",
            "N-arylation",
            "Sonogashira",
            "Heck",
            "Stille",
            "Negishi",
            "Kumada",
            "Ullmann",
        ]

        for rxn_type in coupling_types:
            if checker.check_reaction(rxn_type, rsmi):
                print(f"Detected {rxn_type} coupling reaction: {rsmi}")
                return True

        # Also check for C-C, C-N, C-O bond formation
        reactants_part = rsmi.split(">")[0]
        product_part = rsmi.split(">")[-1]

        # If we have 2 significant reactants and they're combined in the product
        reactants = [r for r in reactants_part.split(".") if not is_reagent_or_catalyst(r)]
        if len(reactants) == 2:
            # This is a simplistic check - in a real implementation, we would do more
            # sophisticated analysis of bond formation between fragments
            return True

        return False

    def dfs_traverse(node):
        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            all_reactants = rsmi.split(">")[0].split(".")

            # Filter out reagents and catalysts
            significant_reactants = [r for r in all_reactants if not is_reagent_or_catalyst(r)]

            # Count number of significant fragments
            fragment_count = len(significant_reactants)
            fragment_counts.append(fragment_count)

            # Check if this is a coupling reaction
            is_coupling = is_coupling_reaction(rsmi)
            coupling_reactions.append(is_coupling)

            print(f"Reaction with {fragment_count} significant fragments: {rsmi}")
            print(f"Is coupling reaction: {is_coupling}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if it's a linear assembly:
    # 1. At least 2 coupling reactions
    # 2. At least 60% of reactions are coupling reactions with 2 fragments
    has_coupling_reactions = len(coupling_reactions) >= 2
    coupling_percentage = (
        coupling_reactions.count(True) / len(coupling_reactions) if coupling_reactions else 0
    )
    two_fragment_percentage = (
        fragment_counts.count(2) / len(fragment_counts) if fragment_counts else 0
    )

    is_linear = (
        has_coupling_reactions and coupling_percentage >= 0.6 and two_fragment_percentage >= 0.6
    )

    print(f"Fragment counts: {fragment_counts}")
    print(f"Coupling reactions: {coupling_reactions}")
    print(
        f"Is linear assembly: {is_linear} (coupling: {coupling_percentage:.2f}, two-fragment: {two_fragment_percentage:.2f})"
    )

    return is_linear
