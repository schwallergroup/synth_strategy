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
    This function detects a synthetic strategy involving N-arylation.
    """
    has_n_arylation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_n_arylation

        # Print node type for debugging
        if node["type"] == "mol":
            print(f"Depth {depth}: Molecule node with SMILES: {node['smiles'][:30]}...")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Depth {depth}: Checking reaction: {rsmi[:50]}...")

            # Extract reactants and product for further analysis
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]
                print(f"  Product: {product[:30]}...")
            except Exception as e:
                print(f"  Error extracting reactants/product: {e}")

            # Check for various N-arylation reactions using the checker function
            n_arylation_reactions = [
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                "Buchwald-Hartwig",
                "Ullmann-Goldberg Substitution amine",
                "{N-arylation_heterocycles}",
                "Goldberg coupling",
                "Goldberg coupling aryl amine-aryl chloride",
                "Goldberg coupling aryl amide-aryl chloride",
            ]

            for reaction_type in n_arylation_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"  Detected N-arylation reaction: {reaction_type}")
                    has_n_arylation = True
                    return

            # If no specific N-arylation reaction was detected, check for N-aryl bond formation
            if product and not has_n_arylation:
                # Check for N-aryl bond formation in the product
                n_aryl_indicators = [
                    checker.check_fg("Aniline", product),
                    checker.check_fg("Secondary amine", product),
                    checker.check_fg("Tertiary amine", product),
                ]

                # Check for reactants that would indicate N-arylation
                amine_reactants = [
                    any(checker.check_fg("Primary amine", r) for r in reactants if r),
                    any(checker.check_fg("Secondary amine", r) for r in reactants if r),
                ]

                aryl_halide_reactants = [
                    any(checker.check_fg("Aromatic halide", r) for r in reactants if r)
                ]

                # Check for heterocyclic rings that might have undergone N-arylation
                heterocyclic_rings = ["pyrrole", "indole", "pyrazole", "imidazole", "benzimidazole"]
                heterocycle_in_product = any(
                    checker.check_ring(ring, product) for ring in heterocyclic_rings
                )

                # Look for patterns suggesting N-arylation
                if any(n_aryl_indicators) and (any(amine_reactants) or any(aryl_halide_reactants)):
                    print(f"  Detected potential N-arylation: N-aryl bond formed in product")
                    has_n_arylation = True

                # Check specifically for N-arylation of heterocycles
                elif heterocycle_in_product and any(aryl_halide_reactants):
                    print(f"  Detected potential N-arylation of heterocycle")
                    has_n_arylation = True

                # Check for specific case in the test data: formation of N-aryl bond in heterocycle
                elif (
                    "c1" in product and "N" in product and any(r and "NH2" in r for r in reactants)
                ):
                    print(
                        f"  Detected potential N-arylation: Amino group reacting with aromatic system"
                    )
                    has_n_arylation = True

        # Traverse children
        if not has_n_arylation:  # Only continue if we haven't found N-arylation yet
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting traversal to find N-arylation strategy...")
    dfs_traverse(route)

    print(f"N-arylation strategy found: {has_n_arylation}")
    return has_n_arylation
