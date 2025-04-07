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
    This function detects late-stage N-alkylation, particularly at the final step of synthesis.
    """
    final_step_is_n_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_n_alkylation

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Check if this is a late-stage step (depth 0 or 1)
            if depth <= 1:
                print(f"Checking for N-alkylation at late stage (depth {depth})")

                # Direct check for N-alkylation reaction types
                if (
                    checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction("Methylation with MeI_primary", rsmi)
                    or checker.check_reaction("Methylation with MeI_secondary", rsmi)
                    or checker.check_reaction("Methylation with MeI_tertiary", rsmi)
                    or checker.check_reaction("Eschweiler-Clarke Primary Amine Methylation", rsmi)
                    or checker.check_reaction("Eschweiler-Clarke Secondary Amine Methylation", rsmi)
                    or checker.check_reaction(
                        "Reductive methylation of primary amine with formaldehyde", rsmi
                    )
                    or checker.check_reaction("N-methylation", rsmi)
                ):

                    print(f"Detected late-stage N-alkylation reaction: {rsmi}")
                    final_step_is_n_alkylation = True
                else:
                    # Fallback check for N-alkylation by examining reactants and products
                    has_primary_amine = False
                    has_secondary_amine = False
                    has_halide = False

                    # Check if product contains secondary or tertiary amine
                    has_secondary_amine_product = checker.check_fg("Secondary amine", product)
                    has_tertiary_amine_product = checker.check_fg("Tertiary amine", product)

                    # Check reactants for amines and halides
                    for reactant in reactants:
                        if checker.check_fg("Primary amine", reactant):
                            has_primary_amine = True
                        if checker.check_fg("Secondary amine", reactant):
                            has_secondary_amine = True

                        if (
                            checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                            or checker.check_fg("Aromatic halide", reactant)
                        ):
                            has_halide = True

                    # Check for appropriate amine transformation
                    if has_primary_amine and has_halide and has_secondary_amine_product:
                        print(f"Detected late-stage N-alkylation of primary amine: {rsmi}")
                        final_step_is_n_alkylation = True
                    elif has_secondary_amine and has_halide and has_tertiary_amine_product:
                        print(f"Detected late-stage N-alkylation of secondary amine: {rsmi}")
                        final_step_is_n_alkylation = True

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    return final_step_is_n_alkylation
