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
    This function detects if the synthetic route involves late-stage functionalization,
    particularly of common functional groups like aniline, phenol, or alcohols in the final steps.
    """
    late_stage_functionalization = False

    # List of functional groups commonly involved in late-stage functionalization
    target_fgs = [
        "Aniline",
        "Phenol",
        "Primary alcohol",
        "Secondary alcohol",
        "Primary amine",
        "Secondary amine",
        "Carboxylic acid",
    ]

    # List of reaction types that are commonly used for late-stage functionalization
    lsf_reactions = [
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of primary amines",
        "Acylation of secondary amines",
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Sulfonamide synthesis (Schotten-Baumann) primary amine",
        "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
        "Esterification of Carboxylic Acids",
        "Williamson Ether Synthesis",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_functionalization

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            # Consider depth 0 or 1 as late-stage
            if depth <= 1:
                try:
                    rsmi = node["metadata"]["rsmi"]
                    reactants_smiles = rsmi.split(">")[0]
                    product_smiles = rsmi.split(">")[-1]

                    # Check if this is a known late-stage functionalization reaction
                    for reaction_type in lsf_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(
                                f"Late-stage functionalization reaction detected at depth {depth}: {reaction_type}"
                            )
                            print(f"Reaction SMILES: {rsmi}")
                            late_stage_functionalization = True
                            return

                    # Check for functionalization of specific functional groups
                    for fg in target_fgs:
                        # Check if the functional group exists in the reactants
                        if checker.check_fg(fg, reactants_smiles):
                            # For aniline specifically, check if it's being functionalized
                            if fg == "Aniline":
                                # Check for common aniline functionalization reactions
                                if (
                                    checker.check_fg("Secondary amide", product_smiles)
                                    or checker.check_fg(
                                        "Tertiary amide", product_smiles
                                    )
                                    or checker.check_fg("Sulfonamide", product_smiles)
                                ):
                                    print(
                                        f"Late-stage {fg} functionalization detected at depth {depth}"
                                    )
                                    print(f"Reaction SMILES: {rsmi}")
                                    late_stage_functionalization = True
                                    return
                            # For other functional groups, check if they're being modified
                            elif not checker.check_fg(fg, product_smiles):
                                print(
                                    f"Late-stage {fg} modification detected at depth {depth}"
                                )
                                print(f"Reaction SMILES: {rsmi}")
                                late_stage_functionalization = True
                                return
                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_functionalization
