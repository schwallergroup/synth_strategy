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
    Detects N-alkylation with alkyl halide in the synthesis route.
    """
    found_n_alkylation = False

    def dfs_traverse(node):
        nonlocal found_n_alkylation

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]

            # Check for N-alkylation reactions using the checker function
            if checker.check_reaction(
                "N-alkylation of primary amines with alkyl halides", rsmi
            ) or checker.check_reaction(
                "N-alkylation of secondary amines with alkyl halides", rsmi
            ):
                print(f"Found N-alkylation reaction: {rsmi}")
                found_n_alkylation = True

            # If specific reaction check fails, try to verify by examining functional groups
            if not found_n_alkylation:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol and all(reactant_mols):
                    # Check for alkyl halide in reactants
                    has_primary_halide = any(
                        checker.check_fg("Primary halide", r) for r in reactants
                    )
                    has_secondary_halide = any(
                        checker.check_fg("Secondary halide", r) for r in reactants
                    )
                    has_tertiary_halide = any(
                        checker.check_fg("Tertiary halide", r) for r in reactants
                    )

                    # Check for amine in reactants
                    has_primary_amine = any(
                        checker.check_fg("Primary amine", r) for r in reactants
                    )
                    has_secondary_amine = any(
                        checker.check_fg("Secondary amine", r) for r in reactants
                    )

                    # Check for alkylated amine in product
                    has_secondary_amine_product = checker.check_fg(
                        "Secondary amine", product
                    )
                    has_tertiary_amine_product = checker.check_fg(
                        "Tertiary amine", product
                    )

                    # Verify N-alkylation occurred
                    if (
                        has_primary_amine
                        and (
                            has_primary_halide
                            or has_secondary_halide
                            or has_tertiary_halide
                        )
                        and has_secondary_amine_product
                    ):
                        print(f"Found N-alkylation of primary amine: {rsmi}")
                        found_n_alkylation = True
                    elif (
                        has_secondary_amine
                        and (
                            has_primary_halide
                            or has_secondary_halide
                            or has_tertiary_halide
                        )
                        and has_tertiary_amine_product
                    ):
                        print(f"Found N-alkylation of secondary amine: {rsmi}")
                        found_n_alkylation = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return found_n_alkylation
