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
    Detects if the synthetic route incorporates a cycloalkyl group via alkylation.
    """
    cycloalkyl_incorporation_found = False

    # List of cycloalkyl rings to check
    cycloalkyl_rings = [
        "cyclopropane",
        "cyclobutane",
        "cyclopentane",
        "cyclohexane",
        "cycloheptane",
        "cyclooctane",
    ]

    # List of alkylation reaction types
    alkylation_reactions = [
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "S-alkylation of thiols",
        "S-alkylation of thiols (ethyl)",
        "S-alkylation of thiols with alcohols",
        "S-alkylation of thiols with alcohols (ethyl)",
        "Williamson Ether Synthesis",
        "Alkylation of amines",
        "Methylation with MeI_primary",
        "Methylation with MeI_secondary",
        "Methylation with MeI_tertiary",
        "Methylation with MeI_aryl",
        "Methylation with MeI_SH",
        "Methylation",
        "Methylation of OH with DMS",
        "Methylation with DMS",
        "Methylation with DMC",
        "DMS COOH methylation",
        "DMS Amine methylation",
        "N-methylation",
        "S-methylation",
        "O-methylation",
        "C-methylation",
        "Friedel-Crafts alkylation",
        "Friedel-Crafts alkylation with halide",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal cycloalkyl_incorporation_found

        if cycloalkyl_incorporation_found:
            return  # Early return if already found

        print(f"Examining node at depth {depth}: {node.get('type', 'unknown')}")

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Reaction SMILES: {rsmi}")

                # In retrosynthetic direction:
                # - The "product" is the starting material (forward direction)
                # - The "reactants" are the retrosynthetic products
                forward_product = rsmi.split(">")[-1]
                forward_reactants = rsmi.split(">")[0].split(".")

                # Check if this is an alkylation reaction
                is_alkylation = False
                for reaction_type in alkylation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_alkylation = True
                        print(f"Found alkylation reaction: {reaction_type}")
                        break

                # Special case for C-C bond formation reactions with organolithium or similar
                if not is_alkylation and ("[Li]" in rsmi or "Li" in rsmi):
                    print("Checking for organolithium-mediated alkylation")
                    # Look for alkyl halide with cycloalkyl
                    for reactant in forward_reactants:
                        for ring in cycloalkyl_rings:
                            if checker.check_ring(ring, reactant) and any(
                                x in reactant for x in ["I", "Br", "Cl"]
                            ):
                                is_alkylation = True
                                print(f"Found organolithium-mediated alkylation with {ring}")
                                break
                        if is_alkylation:
                            break

                # Special case for C-alkylation reactions not covered by standard types
                if not is_alkylation:
                    for reactant in forward_reactants:
                        # Check for alkyl halide with cycloalkyl
                        for ring in cycloalkyl_rings:
                            if checker.check_ring(ring, reactant) and any(
                                x in reactant for x in ["I", "Br", "Cl"]
                            ):
                                # Check if product has the cycloalkyl ring and a new C-C bond
                                if checker.check_ring(ring, forward_product):
                                    is_alkylation = True
                                    print(f"Found C-C bond formation with {ring}")
                                    break
                        if is_alkylation:
                            break

                if is_alkylation:
                    # Check for cycloalkyl group in forward reactants (alkylating agents)
                    cycloalkyl_reactant = None
                    cycloalkyl_ring_found = None

                    for reactant in forward_reactants:
                        for ring in cycloalkyl_rings:
                            if checker.check_ring(ring, reactant):
                                cycloalkyl_reactant = reactant
                                cycloalkyl_ring_found = ring
                                print(f"Found {ring} in reactant: {reactant}")
                                break
                        if cycloalkyl_reactant:
                            break

                    # Check if the cycloalkyl group is incorporated into the forward product
                    if cycloalkyl_reactant and checker.check_ring(
                        cycloalkyl_ring_found, forward_product
                    ):
                        # Verify this is an actual incorporation (not already present)
                        # Check if any of the other reactants already had the cycloalkyl group
                        other_reactants = [r for r in forward_reactants if r != cycloalkyl_reactant]

                        # If there's at least one other reactant without the cycloalkyl group,
                        # and the product has it, then it was incorporated
                        if other_reactants:
                            for other_reactant in other_reactants:
                                if not checker.check_ring(cycloalkyl_ring_found, other_reactant):
                                    print(
                                        f"Found reactant without {cycloalkyl_ring_found}: {other_reactant}"
                                    )
                                    print(
                                        f"Confirmed {cycloalkyl_ring_found} incorporation via alkylation"
                                    )
                                    cycloalkyl_incorporation_found = True
                                    return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {cycloalkyl_incorporation_found}")

    return cycloalkyl_incorporation_found
