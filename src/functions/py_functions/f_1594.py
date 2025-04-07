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
    Detects a late-stage N-alkylation (in the last two steps of synthesis).
    """
    found_n_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_n_alkylation

        if node["type"] == "reaction" and depth <= 1:  # Only check the last two steps
            try:
                # Extract reaction SMILES
                rsmi = node["metadata"]["rsmi"]

                # Check for N-alkylation reactions using the checker
                n_alkylation_reaction_types = [
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "Methylation with MeI_primary",
                    "Methylation with MeI_secondary",
                    "Methylation with MeI_tertiary",
                    "DMS Amine methylation",
                    "Eschweiler-Clarke Primary Amine Methylation",
                    "Eschweiler-Clarke Secondary Amine Methylation",
                    "Reductive methylation of primary amine with formaldehyde",
                    "N-methylation",
                    "Alkylation of amines",
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "Reductive amination with alcohol",
                ]

                for reaction_type in n_alkylation_reaction_types:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Detected late-stage N-alkylation: {reaction_type} at depth {depth}"
                        )
                        found_n_alkylation = True
                        return  # Exit early once found

                # If no specific reaction type matched, check for general pattern
                if not found_n_alkylation:
                    # Extract reactants and product
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Check if any reactant has primary or secondary amine
                    has_amine_reactant = False
                    for r_smiles in reactants_smiles:
                        if checker.check_fg(
                            "Primary amine", r_smiles
                        ) or checker.check_fg("Secondary amine", r_smiles):
                            has_amine_reactant = True
                            break

                    # Check if any reactant has alkylating agent
                    has_alkylating_agent = False
                    for r_smiles in reactants_smiles:
                        if (
                            checker.check_fg("Primary halide", r_smiles)
                            or checker.check_fg("Secondary halide", r_smiles)
                            or checker.check_fg("Tertiary halide", r_smiles)
                            or checker.check_fg("Mesylate", r_smiles)
                            or checker.check_fg("Tosylate", r_smiles)
                            or checker.check_fg("Triflate", r_smiles)
                        ):
                            has_alkylating_agent = True
                            break

                    # Check for carbonyl compounds for reductive amination
                    has_carbonyl = False
                    for r_smiles in reactants_smiles:
                        if (
                            checker.check_fg("Aldehyde", r_smiles)
                            or checker.check_fg("Ketone", r_smiles)
                            or checker.check_fg("Formaldehyde", r_smiles)
                        ):
                            has_carbonyl = True
                            break

                    # Count tertiary amines in reactants and product
                    reactant_tertiary_count = sum(
                        1
                        for r in reactants_smiles
                        if checker.check_fg("Tertiary amine", r)
                    )
                    product_has_tertiary = checker.check_fg(
                        "Tertiary amine", product_smiles
                    )

                    # If we have amine reactant and alkylating agent or carbonyl compound,
                    # and the product has tertiary amine (indicating N-alkylation occurred),
                    # mark as found
                    if (
                        has_amine_reactant
                        and (has_alkylating_agent or has_carbonyl)
                        and product_has_tertiary
                    ):
                        # Additional check: product should have more tertiary amines than reactants combined
                        # or the reaction should convert primary/secondary to tertiary
                        if reactant_tertiary_count == 0 or has_alkylating_agent:
                            print(
                                f"Detected late-stage N-alkylation through functional group analysis at depth {depth}"
                            )
                            found_n_alkylation = True
                            return  # Exit early once found
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage N-alkylation detected: {found_n_alkylation}")
    return found_n_alkylation
