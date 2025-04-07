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
    """Check for late-stage N-alkylation in the synthesis route"""
    found = False
    total_depth = get_max_depth(route)

    def dfs(node, depth=0):
        nonlocal found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for N-alkylation reactions
            if (
                checker.check_reaction(
                    "N-alkylation of primary amines with alkyl halides", rxn_smiles
                )
                or checker.check_reaction(
                    "N-alkylation of secondary amines with alkyl halides", rxn_smiles
                )
                or checker.check_reaction("Alkylation of amines", rxn_smiles)
                or checker.check_reaction("N-methylation", rxn_smiles)
                or checker.check_reaction("Methylation", rxn_smiles)
                or checker.check_reaction("Methylation with MeI_primary", rxn_smiles)
                or checker.check_reaction("Methylation with MeI_secondary", rxn_smiles)
            ):

                # Late stage means low depth (close to final product)
                # Consider first half of the synthesis as late stage
                late_stage_threshold = max(2, total_depth // 2)
                if depth <= late_stage_threshold:
                    # Verify N-alkylation
                    try:
                        reactants = rxn_smiles.split(">")[0].split(".")
                        product = rxn_smiles.split(">")[-1]

                        has_amine_reactant = any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            or checker.check_fg("Aniline", r)
                            for r in reactants
                        )
                        has_alkyl_halide = any(
                            checker.check_fg("Primary halide", r)
                            or checker.check_fg("Secondary halide", r)
                            or checker.check_fg("Tertiary halide", r)
                            for r in reactants
                        )

                        if has_amine_reactant and (
                            has_alkyl_halide
                            or any("C-[Cl,Br,I,F]" in r for r in reactants)
                        ):
                            found = True
                            print(
                                f"Found late-stage N-alkylation at depth {depth} (threshold: {late_stage_threshold}): {rxn_smiles}"
                            )
                    except Exception as e:
                        print(f"Error checking N-alkylation: {e}")
                        # If we can't verify, still consider it found if the reaction type matches
                        found = True
                        print(
                            f"Found late-stage N-alkylation (unverified) at depth {depth}: {rxn_smiles}"
                        )

        # Recursively check children
        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)
    return found
